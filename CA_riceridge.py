import numpy as np
from matplotlib import pyplot as plt
import rasterio as rio
from matplotlib import colors
from itertools import chain
import multiprocessing as mp
import collections
import sys
import cProfile
import subprocess
import netCDF4 as nc
import os



class fireCA:
    
    def __init__(self, path_to_data=False):
        self.plotting = False
        self.show_plot = False
        self.rng = np.random.default_rng()
        
        # constants
        
        self.consts = {'spatial_resolution':250,
                       's_base':2**0,
                       'weight':0.32,
                       'weight_d':1.6,
                       'weight_phi':0.7,
                       'spread_threshold':0.925,
                       'burn_speed':0.07,
                       'spot_sigma':1,
                       'brand_sigma':1,
                       'shrub_speed':1,
                       'forest_speed':1,
                       'smoldering_threshold':20,
                       'extinguished_threshold':5,
                       'ign_loc':[(93,50)],
                       'ign_prob':False,
                       'fname':'',
                       'start_time':np.datetime64('2017-07-24','h')}
                       
        # if it gets pointed to data, initialize with data
        if path_to_data:
            self.get_data(path_to_data)
            self.preprocess()
    
    def get_data(self, path_to_data):
        # ignition probability dataset
        #self.prob_cube = rio.open(path_to_data+'prob.tif').read()
        
        # rate of spread dataset
        self.RoS_cube = np.power(10,rio.open(path_to_data+'RoS.tif').read())
        
        # topography
        self.slope = rio.open(path_to_data+'slope.tif').read()[0]
        self.aspect = rio.open(path_to_data+'aspect.tif').read()[0]
        self.dem = rio.open(path_to_data+'dem.tif').read()[0]
        self.hillshade = rio.open(path_to_data+'hillshade.tif').read()[0]
        
        # wind
        self.wfrom = rio.open(path_to_data+'wdir.tif').read()
        self.wspd = rio.open(path_to_data+'wspd.tif').read()
        
        # vegetation
        self.veg_mtx = rio.open(path_to_data+'veg.tif').read()/100
        self.noveg_mtx = rio.open(path_to_data+'noveg.tif').read()[0]
        
        # initial state
        self.initial_state = (rio.open(path_to_data+'state.tif').read()[0]).astype(int)
        
        # burnday
        self.burnday = np.rint(rio.open(path_to_data+'burnday.tif').read()[0])
    
    def preprocess(self):
        t,self.m,self.n = self.RoS_cube.shape
        self.t = min(self.wfrom.shape[0], t*24)
        
        self.LEFT = 0# int(3/50*self.n)
        self.RIGHT = self.n# int(1/3*self.n)
        self.TOP = 0# int(1/3*self.m)
        self.BOTTOM = self.m# int(2/3*self.m)
        
        self.prob_cube = np.ones_like(self.RoS_cube)*self.consts['ign_prob']
        
        self.burnday[self.burnday<=205] = np.nan
        burn_day = self.burnday[self.TOP:self.BOTTOM,self.LEFT:self.RIGHT]
        burn_day = burn_day[~np.isnan(burn_day)]
        self.VIIRS_burnt_pixels = np.zeros(self.t//24)
        for i in range(self.t//24):
            j = i + 205
            self.VIIRS_burnt_pixels[i] = np.sum(burn_day <= j)
        
        self.upslope = ((270 - self.aspect)*np.pi/180 + np.pi) % (2*np.pi) - np.pi
        
        self.wfrom[self.wfrom<0] = np.nan
        self.wdir = ((270 - self.wfrom)*np.pi/180 + np.pi) % (2*np.pi) - np.pi
        self.wspd[self.wspd<0] = 0
        
        # maximum rate and direction of spread
        # some bullshit i made up, just graph it
        # based off some article i read about the 10% rule
        wm = (self.wspd*3.6) # m s-1 -> km hr-1
        wm = np.exp(0.05*wm)/(1+np.exp(0.6*(10-wm)))
        # slope from Butler, Anderson, Catchpole 2007
        sm = 0.001*np.exp(0.38*self.slope) + 1.6
        # turn area into length
        # area = RoS*4046/62500
        # a = (1-np.exp(-2*np.sqrt(3)*np.pi))/2/np.sqrt(3)
        # b = 2*np.sqrt(2)*(1-np.exp(-np.sqrt(3)*np.pi))/np.sqrt(3)
        # k = 250*(np.sqrt(b**2+4*a*area) - b)/2/a
        k = self.consts['spatial_resolution']*1.732083333553463*(np.sqrt(2.6436051667022675 + 1.1546788547967213*self.RoS_cube*4046/self.consts['spatial_resolution']**2/24) - 1.6259167157952057)
        self.r_max_cube = np.sqrt(np.square(wm) + np.square(sm) + 2*wm*sm*np.cos(self.wdir-self.upslope))*np.repeat(k,24,axis=0)[:wm.shape[0]]
        self.phi_max_cube = (self.upslope + np.arctan2(wm*np.sin(self.wdir-self.upslope), sm + wm*np.cos(self.wdir-self.upslope)))
        
        self.fuel_burn_speed = self.consts['shrub_speed']*self.veg_mtx[0]+self.consts['forest_speed']*self.veg_mtx[1]
        
        self.consts['ign_loc'] = set([tuple(x) for x in np.argwhere(self.initial_state==2)])
        
        self.reset()
    
    def reset(self):
        self.burning = set(self.consts['ign_loc'])
        self.smoldering = set()
        self.extinguished = set()
        self.burnt = set()
        self.unburnable = set([tuple(x) for x in np.argwhere(self.noveg_mtx)])
        self.ign_dict = collections.defaultdict(float)
        self.surface_fuel_dict = collections.defaultdict(float)
        self.bulk_fuel_dict = collections.defaultdict(float)
        
        self.CA_burnt_pixels = np.zeros(self.t-1)
        self.n_prev_burnt = 0
        self.i_prev = 0
        self.s_per_frame = []
        
        self.phi_max = self.phi_max_cube[0]
        self.r_max = self.r_max_cube[0]
        self.wind_spd = self.wspd[0]
        self.wind_dir = self.wdir[0]
        self.prob = self.prob_cube[0]
        self.RoS = self.RoS_cube[0]
    
    
    def plot(self, path=False):
        state_mtx = np.ones_like(self.RoS).astype(int)*3
        for pixel in self.smoldering:
            state_mtx[pixel] = 1
        for pixel in self.burning:
            state_mtx[pixel] = 2
        for pixel in self.burnt:
            state_mtx[pixel] = -1
        for pixel in self.unburnable:
            state_mtx[pixel] = 0
        
        relevant_pixels = list(chain(self.burning, self.smoldering, self.extinguished, self.burnt))
        y = [pixel[0] for pixel in relevant_pixels]
        x = [pixel[1] for pixel in relevant_pixels]
        top = np.min(y)
        bottom = np.max(y)
        left = np.min(x)
        right = np.max(x)
        xpad = max((right-left)/3, 20)
        ypad = max((bottom-top)/3, 20)
        top = int(max(0, top-ypad))
        bottom = int(min(self.m, bottom+ypad))
        left = int(max(0, left-xpad))
        right = int(min(self.n, right+xpad))
        
        #cmap = colors.ListedColormap(['black', 'dimgray', 'red', 'gold', 'darkgreen'])
        cmap = colors.ListedColormap(['black', (0,0,0,0), 'red', 'gold', (0,0,0,0)])
        norm = colors.BoundaryNorm([-1,0,1,2,3,4], cmap.N)
        fig,ax = plt.subplots(1,2,figsize=(15,7))
        ax[0].imshow(self.hillshade[top:bottom,left:right], alpha=1, cmap='gray', interpolation='none', zorder=0)
        ax[0].imshow(state_mtx[top:bottom,left:right], cmap=cmap, alpha=0.6, norm=norm, interpolation='none', zorder=10)
        viirs = np.argwhere(self.burnday[top:bottom,left:right]<=self.day+205).T
        ax[0].scatter(viirs[1], viirs[0], c='r', s=1, alpha=0.4, zorder=5)
        ax[0].tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelright=False, labeltop=False, labelbottom=False)
        ax[1].plot(self.CA_burnt_pixels[:self.day], label='CA')
        ax[1].plot(self.VIIRS_burnt_pixels[:self.day], label='VIIRS')
        ax[1].legend()
        fig.suptitle(f"{str((self.consts['start_time']+self.hour+24*self.day)).replace('T',' ')}:00\n"\
                  f"area burned: {((self.consts['spatial_resolution']/1000)**2*(len(self.burning)+len(self.smoldering)+len(self.burnt))):.2f} km2")
        if path:
            plt.savefig(path)
        if self.show_plot:
            plt.show()
        plt.close()
    
    def step_size(self, s, max_spread):
        if s < 2:
            return 1
        if max_spread*s/self.consts['s_base'] > self.consts['spatial_resolution']:
            s = self.step_size(s/2, max_spread)
        return int(s)
    
    def get_adj_dw(self, pixel):
        weights = {}
        i,j = pixel
        
        r_max = self.r_max[pixel]#*self.s/self.consts['s_base']
        # this line makes the radius nonlinear on RoS. rough substitution for spotting behaviour.
        #r = (r_max/self.consts['spatial_resolution'])**(1 + 1/(1+np.exp(-r_max+100)) + 1/(1+np.exp(-r_max+250))) # in units of grid cells
        r = r_max/self.consts['spatial_resolution'] + 1
        phi = self.phi_max[pixel]
        l = int(r)
        if l < 1: l = 0
        candidates = set([(y,x) for y in range(max(i-1-l, 0),min(i+2+l, self.m)) for x in range(max(j-1-l, 0),min(j+2+l, self.n))]) - set([pixel])
        
        for neighbor in candidates:
            # theta is the direction from the current pixel to the neighbor
            theta = np.arctan2(pixel[0]-neighbor[0],neighbor[1]-pixel[1])
            diff = np.abs(phi - theta)
            d = np.linalg.norm(np.subtract(pixel,neighbor))
            
            # the first condition defines the neighborhood shape
            if (d <= np.sqrt(2)+r*np.exp(-np.sqrt(3)*diff)) and (~np.isnan(self.r_max[neighbor]))\
            and (neighbor not in chain(self.burning,self.unignitable,self.burnt,self.unburnable)):
                # weights are 2-d gaussian in distance and angle, scaled by RoS
                weights[neighbor] = self.r_max[neighbor]/100*self.consts['weight']*np.exp((-(d/((r+1)*self.consts['weight_d']))**2-(diff/(np.pi/2*self.consts['weight_phi']))**2)/2)
        return weights
    
    def smolder(self):
        # progress burn in all smoldering cells
        # if RoS falls below threshold, transition to extinguished
        # if any smoldering cells reach the threshold, transition them to extinguished
        remove_from_smoldering = set()
        for pixel in self.smoldering:
            if self.RoS[pixel] <= self.consts['extinguished_threshold']:
                remove_from_smoldering.add(pixel)
                self.extinguished.add(pixel)
            elif (self.RoS[pixel] > self.consts['smoldering_threshold'] and self.surface_fuel_dict[pixel]):
                remove_from_smoldering.add(pixel)
                self.burning.add(pixel)
            else:
                spread = self.RoS[pixel]*self.s/self.consts['s_base']*self.consts['burn_speed']
                self.bulk_fuel_dict[pixel] += spread
                if self.bulk_fuel_dict[pixel] >= self.consts['spatial_resolution']:
                    del self.bulk_fuel_dict[pixel]
                    remove_from_smoldering.add(pixel)
                    self.burnt.add(pixel)
        self.smoldering.difference_update(remove_from_smoldering)
    
    def spread(self):
        # progress fire in all burning cells
        # if RoS falls below threshold, transition to smoldering
        # if any burning cells reach the threshold, transition them to smoldering
        # along the way, release firebrands
        remove_from_burning = set()
        for pixel in self.burning:
            if self.RoS[pixel] <= self.consts['extinguished_threshold']:
                remove_from_burning.add(pixel)
                self.extinguished.add(pixel)
            elif self.RoS[pixel] <= self.consts['smoldering_threshold']:
                remove_from_burning.add(pixel)
                self.smoldering.add(pixel)
            else:
                spread = self.r_max[pixel]*self.s/self.consts['s_base']*self.consts['burn_speed']*self.fuel_burn_speed[pixel]
                self.surface_fuel_dict[pixel] += spread
                if self.surface_fuel_dict[pixel] >= self.consts['spatial_resolution']:
                    del self.surface_fuel_dict[pixel]
                    remove_from_burning.add(pixel)
                    self.smoldering.add(pixel)
                """
                ignited = set()
                # forest fraction times meters spread somehow becomes firebrands released
                n_brands = int(self.veg_mtx[5][pixel]*spread)
                # each firebrand samples from a distribution and that number is the distance it travels
                # the distribution should depend on wind speed and intensity (sub RoS)
                sigma = (self.wind_spd[pixel]*2.24/50*1609/self.consts['spatial_resolution'])*self.consts['spot_sigma']
                reach = np.abs(self.rng.normal(0,sigma,n_brands))
                reach = reach[reach>1]
                phi = self.wind_dir[pixel]
                brands = [(y,x) for x in np.rint(pixel[1]+reach*np.cos(phi)).astype(int) for y in np.rint(pixel[0]-reach*np.sin(phi)).astype(int)]
                branded = set(brands)
                # if the firebrand lands in an ignitable cell, that cell ignites (with a probability?)
                for neighbor in branded:
                    if (0 <= neighbor[0] < self.m) and (0 <= neighbor[1] < self.n) and (~np.isnan(self.r_max[neighbor]))\
                    and (neighbor not in chain(self.burning,self.smoldering,self.burnt,self.unburnable))\
                    and (np.any(self.rng.normal(0,self.consts['brand_sigma'],np.sum([1 for coord in brands if coord == neighbor]))>1-self.prob[neighbor])):
                        ignited.update(neighbor)
                self.burning.update(ignited)
                """
        self.burning.difference_update(remove_from_burning)
    
    def ignite(self):
        spreading = set()
        for pixel in self.burning:
            if self.surface_fuel_dict[pixel] >= self.consts['spatial_resolution']*self.consts['spread_threshold']:
                i,j = pixel
                neighbors = set([(y,x) for y in range(max(i-1, 0),min(i+2, self.m)) for x in range(max(j-1, 0),min(j+2, self.n))]) - set([pixel])
                #if len(neighbors.intersection(self.burning)) != len(neighbors):
                if True:
                    spreading.update([pixel])
        
        self.unignitable = set([pixel for pixel in chain(self.extinguished,self.smoldering) if self.surface_fuel_dict[pixel]-0<=1e-12])
        if len(spreading) > 10:
            with mp.Pool(16) as p:
                weight_list = [x.items() for x in p.map(self.get_adj_dw, spreading)]
                p.close()
                p.join()
        else:
            weight_list = [x.items() for x in [self.get_adj_dw(pixel) for pixel in spreading]]
        
        try:
            for neighbor,weight in chain(*weight_list):
                self.ign_dict[neighbor] += weight * self.s/self.consts['s_base']
                if self.ign_dict[neighbor] >= 1-self.prob[neighbor]:
                    del self.ign_dict[neighbor]
                    self.burning.update([neighbor])
        except ValueError:
            pass
    
    def run(self):
        # take t and multiply it by 32 (s_base)
        # keep track of how many 32ths have passed
        # begin a while loop
        # update weather conditions
        # determine step size
        # transition smoldering to burnt
        # progress fire in all burning cells
        # if any burning cells reach the threshold, transition them to smoldering
        # ignite susceptible cells
        # increment 32ths
        
        self.i = 0
        self.hour = 0
        self.day = 0
        while self.i < (self.t-1)*self.consts['s_base']:
            
            # exit condition
            if len(self.burning) > 1000:
                print("Number of burning cells exceeds 1000")
                self.plot(self.consts['gif_dir']+f"{self.day:02}_{self.hour:02}:{self.minute:02}:{self.second:02}.png")
                break
                
            
            self.n_burnt = len(self.burning)+len(self.burnt)+len(self.smoldering)+len(self.extinguished)
            
            # update weather conditions
            if self.i > self.consts['s_base']*(24*self.day+self.hour):
                j = int(self.i//self.consts['s_base'])
                if self.hour != j%24:
                    self.hour = j%24
                    self.wind_spd = self.wspd[j]
                    self.wind_dir = self.wdir[j]
                    self.r_max = self.r_max_cube[j]
                    self.phi_max = self.phi_max_cube[j]
                    if self.hour == 0:
                        self.day = int(j//24)
                        self.CA_burnt_pixels[self.day] = self.n_burnt
                        self.RoS = self.RoS_cube[self.day]
                        self.prob = self.prob_cube[self.day]
            
            if ((self.plotting) and ((self.n_burnt - self.n_prev_burnt > 0.05*self.n_prev_burnt) or (self.i - self.i_prev > 6*self.consts['s_base']))):
                self.n_prev_burnt = self.n_burnt
                self.i_prev = self.i
                self.minute = int(((self.i/self.consts['s_base'])%1)*60)
                self.second = int(((((self.i/self.consts['s_base'])%1)*60)%1)*60)
                self.plot(self.consts['gif_dir']+f"{self.day:02}_{self.hour:02}:{self.minute:02}:{self.second:02}.png")
                self.s_per_frame.append(self.i - self.i_prev)
            
            # determine step size
            try:
                max_spread = np.max([self.r_max[pixel] for pixel in self.burning])
            except ValueError:
                max_spread = 0
                if len(self.smoldering) == 0:
                    print(f'Fire extinguished.')
                    self.plot(self.consts['gif_dir']+f"{self.day:02}_{self.hour:02}:{self.minute:02}:{self.second:02}.png")
                    break
            self.s = self.step_size(self.consts['s_base'], max_spread)
            
            # smolder
            self.smolder()
            
            # progress fire in all burning cells
            self.spread()
            
            # ignite neighboring cells
            self.ignite()
            
            # heat escapes from cells which are heating up
            self.ign_dict.update((x,y*(1-self.s/self.consts['s_base'])) for x,y in self.ign_dict.items())
            
            # increment 32ths
            self.i += self.s
            
def main():
    a = sys.argv[1:]
    print("initializing CA")
    CA = fireCA(a[0])
    CA.plotting = True
    CA.consts['gif_dir'] = a[1]
    for key,val in np.loadtxt(a[2], delimiter=',',dtype=str).T:
        CA.consts[key] = float(val)
    CA.consts['spatial_resolution'] = int(a[3])
    #CA.consts['weight_d'] *= (240/CA.consts['spatial_resolution'])
    CA.consts['burn_speed'] *= (CA.consts['spatial_resolution']/240)**np.sqrt(3)
    with cProfile.Profile() as p:
        print("running CA")
        p.runcall(CA.run)
        p.dump_stats(f'CAprofile_{a[4]}')
    np.savetxt(a[1]+'burnt.csv', CA.CA_burnt_pixels, delimiter=',')
    
    print("making gif")
    file_list = []
    ffconcat = ["ffconcat version 1.0\n"]
    directory = '/home/adam.viray/Documents/CA/output/riceridge/' + a[3] + 'm/'
    os.chdir(directory)
    for root,dirs,files in os.walk(directory):
        for filename in files:
            if filename.endswith(".png") and len(filename) == 15:
                file_list.append("file " + filename + "\n")
    file_list.sort()
    datetime_list = [np.datetime64(str(np.datetime64(int(fname[5:-5].split("_")[0]),'D'))+"T"+fname[5:-5].split("_")[1]) for fname in file_list]
    duration_list = np.append(np.diff(datetime_list).astype('timedelta64[s]'), [3600]).astype(float)*32/3600/240
    duration_list = [f"duration {duration:.3f}\n" for duration in duration_list]
    ffconcat += [item for pair in zip(file_list, duration_list) for item in pair]
    ffconcat += ffconcat[-2:]
    fffilename = 'ffconcat.txt'
    with open(fffilename, 'w') as f:
        f.writelines(ffconcat)
    subprocess.run(f"ffmpeg -y -f concat -safe 0 -i {fffilename} {a[3]}m.gif".split(" "))
    #subprocess.run(f"ffmpeg -y -f concat -safe 0 -i {fffilename} -c:v libx264 -pix_fmt yuv420p -movflags +faststart output.mp4".split(" "))

if __name__=='__main__':
    main()
