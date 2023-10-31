import numpy as np
from itertools import chain
import sys
import rasterio as rio
import matplotlib.pyplot as plt
from matplotlib import colors

class fireCA:
    
    def __init__(self, path_to_data=False):
        self.plotting = False
        self.rng = np.random.default_rng()
        
        # constants
        self.consts = {'spatial_resolution':250,
                       'veg_dens_threshold':0.25,
                       'bare_dirt_threshold':0.1,
                       's_base':2**7,
                       'weight':0.325,
                       'weight_d':1.6,
                       'weight_phi':0.7,
                       'spread_threshold':0.925,
                       'burn_speed':0.07,
                       'spot_sigma':0.3,
                       'brand_sigma':0.3,
                       'fname':''}
                       
        # if it gets pointed to data, initialize with data
        if path_to_data:
            self.get_data(path_to_data)
            self.preprocess()
        
    def get_data(self, path_to_data):
        
        # ignition probability dataset
        self.prob_cube = rio.open(path_to_data+'prob.tif').read()
        
        # rate of spread dataset
        self.RoS_cube = rio.open(path_to_data+'RoS.tif').read()
        
        self.slope = rio.open(path_to_data+'slope.tif').read()[0]
        self.aspect = rio.open(path_to_data+'aspect.tif').read()[0]
        self.dem = rio.open(path_to_data+'dem.tif').read()[0]
        
        # wind
        self.wfrom = rio.open(path_to_data+'wdir.tif').read()
        self.wspd = rio.open(path_to_data+'wspd.tif').read()
        
        # vegetation
        self.veg_mtx = rio.open(path_to_data+'veg.tif').read()/100
        
        # initial state
        self.initial_state = rio.open(path_to_data+'state.tif').read()[0]
        
        
    def preprocess(self):
        t,self.m,self.n = self.prob_cube.shape
        self.t = min(self.wfrom.shape[0], t*24)
        
        self.upslope = ((270 - self.aspect)*np.pi/180 + np.pi) % (2*np.pi) - np.pi
        
        self.wfrom[self.wfrom<0] = np.nan
        self.wdir = ((270-self.wfrom)*np.pi/180 + np.pi) % (2*np.pi) - np.pi
        self.wspd[self.wspd<0] = 0
        
        # maximum rate and direction of spread
        # some bullshit i made up, just graph it
        wm = (self.wspd*3.6)
        wm = np.exp(0.05*wm)/(1+np.exp(0.6*(10-wm)))
        # slope from Butler, Anderson, Catchpole 2007
        sm = 0.001*np.exp(0.38*np.arctan(self.slope/100)*180/np.pi) + 1.6
        # turn area into length
        # area = RoS*4046/62500
        # a = (1-np.exp(-2*np.sqrt(3)*np.pi))/2/np.sqrt(3)
        # b = 2*np.sqrt(2)*(1-np.exp(-np.sqrt(3)*np.pi))/np.sqrt(3)
        # k = 250*(np.sqrt(b**2+4*a*area) - b)/2/a
        k = 433.02083338836576*(np.sqrt(2.6436051667022675 + 1.1546788547967213*self.RoS_cube*4046/62500/24) - 1.6259167157952057)
        self.r_max_cube = np.sqrt(np.square(wm) + np.square(sm) + 2*wm*sm*np.cos(self.wdir-self.upslope))*np.repeat(k,24,axis=0)[:wm.shape[0]]
        self.phi_max_cube = (self.upslope + np.arctan2(wm*np.sin(self.wdir-self.upslope), sm + wm*np.cos(self.wdir-self.upslope)))
        
        self.fuel_burn_speed = 0.8*self.veg_mtx[0]+1.2*self.veg_mtx[2]
        
        self.consts['ign_loc'] = set([tuple(x) for x in np.argwhere(self.initial_state==2)])
        
        self.reset()
        
    def reset(self):
        self.burning = set(self.consts['ign_loc'])
        self.smoldering = set()
        self.burnt = set()
        self.unburnable = set([tuple(x) for x in np.argwhere(self.veg_mtx[0]+self.veg_mtx[2]<self.consts['veg_dens_threshold'])])\
                          .union(set([tuple(x) for x in np.argwhere(self.veg_mtx[1]>self.consts['bare_dirt_threshold'])]))
        self.ign_dict = {}
        self.burn_dict = {}
        
        self.CA_burnt_pixels = np.zeros(self.t-1)
        
        self.phi_max = self.phi_max_cube[0]
        self.r_max = self.r_max_cube[0]
        self.wind_spd = self.wspd[0]
        self.wind_dir = self.wdir[0]
        self.RoS = self.RoS_cube[0]
        self.prob = self.prob_cube[0]
        
        
    def plot(self, path=False, show=True):
        state_mtx = np.ones_like(self.prob).astype(int)*3
        for pixel in self.smoldering:
            state_mtx[pixel] = 1
        for pixel in self.burning:
            state_mtx[pixel] = 2
        for pixel in self.burnt:
            state_mtx[pixel] = -1
        for pixel in self.unburnable:
            state_mtx[pixel] = 0
        
        cmap = colors.ListedColormap(['black', 'dimgray', 'red', 'gold', 'darkgreen'])
        norm = colors.BoundaryNorm([-1,0,1,2,3,4], cmap.N)
        plt.imshow(state_mtx, cmap=cmap, norm=norm, interpolation='none', zorder=10)
        plt.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelright=False, labeltop=False, labelbottom=False)
        #plt.imshow(self.dem, alpha=0.4, interpolation='none', zorder=0)
        plt.title(f'day {self.day} hour {self.hour}')
        #plt.title(f"""
        #             veg_dens_threshold: {self.consts['veg_dens_threshold']}, bare_dirt_threshold: {self.consts['bare_dirt_threshold']},
        #             weight: {self.consts['weight']}, weight_d: {self.consts['weight_d']}, weight_phi: {self.consts['weight_phi']},
        #             burn_speed:{self.consts['burn_speed']}, spread_threshold: {self.consts['spread_threshold']},
        #             \n\n day {self.day} hour {self.hour}""", y=1)
        if path:
            plt.savefig(path)
        if show:
            plt.show()
    
    def step(self, s, max_spread):
        if s < 2:
            return 1
        if max_spread*s/self.consts['s_base'] > self.consts['spatial_resolution']:
            s = self.step(s/2, max_spread)
        return int(s)
    
    def get_adj_dw(self, pixel):
        neighbors = set()
        weights = {}
        i,j = pixel
        
        r_max = self.r_max[pixel]
        # this line makes the radius nonlinear on RoS. rough substitution for spotting behaviour.
        #r = (r_max/self.consts['spatial_resolution'])**(1 + 1/(1+np.exp(-r_max+100)) + 1/(1+np.exp(-r_max+250))) # in units of grid cells
        r = r_max/self.consts['spatial_resolution'] + 1
        phi = self.phi_max[pixel]
        
        l = int(r)
        if l < 1: l = 0
        candidates = set([(y,x) for y in range(i-1-l,i+2+l) for x in range(j-1-l,j+2+l)]) - set([pixel])

        for neighbor in candidates:
            # theta is the direction from the current pixel to the neighbor
            theta = np.arctan2(pixel[0]-neighbor[0],neighbor[1]-pixel[1])
            diff = np.abs(phi - theta)
            d = np.linalg.norm(np.subtract(pixel,neighbor))
            
            # the first line defines the neighborhood shape
            if (d <= np.sqrt(2)+r*np.exp(-np.sqrt(3)*diff))\
            and (0 <= neighbor[0] < self.m) and (0 <= neighbor[1] < self.n) and (~np.isnan(self.r_max[neighbor]))\
            and (neighbor not in chain(self.burning,self.smoldering,self.burnt,self.unburnable)):
                # weights are 2-d gaussian in distance and angle, scaled by RoS
                weights[neighbor] = self.r_max[neighbor]/50*self.consts['weight']*np.exp((-(d/((r+1)*self.consts['weight_d']))**2-(diff/(np.pi/2*self.consts['weight_phi']))**2)/2)
                neighbors.add(neighbor)
        
        return neighbors, weights
    
    def burn(self, pixel, p):
        if pixel not in self.ign_dict.keys():
            self.ign_dict[pixel] = 0
        self.ign_dict[pixel] += p
        if self.ign_dict[pixel] >= 1-self.prob[pixel]:
            del self.ign_dict[pixel]
            return [pixel]
        return []
    
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
            
            # update weather conditions
            if self.i > self.consts['s_base']*(24*self.day+self.hour):
                j = int(self.i//self.consts['s_base'])
                if self.hour != j%24:
                    if self.plotting:
                        self.plot(self.consts['gif_dir']+f'{self.day:02}_{self.hour:02}.png',False)
                    self.hour = j%24
                    self.wind_spd = self.wspd[j]
                    self.wind_dir = self.wdir[j]
                    self.r_max = self.r_max_cube[j]
                    self.phi_max = self.phi_max_cube[j]
                    if self.hour == 0:
                        self.day = int(j//24)
                        self.prob = self.prob_cube[self.day]
                        self.CA_burnt_pixels[self.day] = len(self.burning)+len(self.smoldering)+len(self.burnt)
                        


            # determine step size
            try:
                max_spread = np.max([self.r_max[pixel] for pixel in self.burning])
            except ValueError:
                max_spread = 0
                if len(self.smoldering) == 0:
                    print(f'Fire extinguished.')
                    break
            self.s = self.step(self.consts['s_base'], max_spread)
            
            # heat escapes from cells which are heating up
            self.ign_dict.update((x,y*(1-self.s/self.consts['s_base'])) for x,y in self.ign_dict.items())
            
            # transition smoldering to burnt
            self.burnt.update(self.smoldering)
            self.smoldering = set()
            
            
            ignited = set()
            
            # progress fire in all burning cells
            # if any burning cells reach the threshold, transition them to smoldering
            # along the way, release firebrands
            remove_from_burning = set()
            for pixel in self.burning:
                if pixel not in self.burn_dict.keys():
                    self.burn_dict[pixel] = 0
                spread = self.r_max[pixel]*self.s/self.consts['s_base']*self.consts['burn_speed']*self.fuel_burn_speed[pixel]
                self.burn_dict[pixel] += spread
                if self.burn_dict[pixel] >= self.consts['spatial_resolution']:
                    remove_from_burning.add(pixel)
                    self.smoldering.add(pixel)
                """
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
                """
            self.burning.difference_update(remove_from_burning)

            # ignite neighboring cells
            for pixel in self.burning:
                if self.burn_dict[pixel] >= self.consts['spatial_resolution']*self.consts['spread_threshold']:
                    adj,weight = self.get_adj_dw(pixel)
                    for neighbor in adj:
                        p = weight[neighbor] * self.s/self.consts['s_base']
                        ignited.update(self.burn(neighbor, p))

            """
            for pixel in chain(self.smoldering):
                adj,weight = self.get_adj_dw(pixel)
                for neighbor in adj:
                    p = weight[neighbor] * self.s/self.consts['s_base']
                    ignited.update(self.burn(neighbor, p))
            """
            self.burning.update(ignited)


            # increment 32ths
            self.i += self.s

# arguments should be:
# 1. path to VIIRS folder. should look like /mnt/mordor3/data/VIIRS/riceridge/model_inputs/
# 2. path to output file. should look like /home/adam.viray/Documents/CA/output/1.csv
# 3. parameter 'weight'. default 0.5, could be between 0 and 10
# 4. parameter 'weight_d'. default 0.8, could be between 0 and 10
# 5. parameter 'weight_phi'. default 0.8, could be between 0 and 1
# 6. parameter 'spread_threshold'. default 0.95, could be between 0 and 1
# 7. parameter 'burn_speed'. default 0.06, could be between 0 and 10
# 8. parameter 'spot_sigma'.
# 9. parameter 'brand_sigma'.

def main():
    a = sys.argv[1:]
    CA = fireCA(a[0])
    
    if len(a) > 2:
        CA.consts.update({
            'weight':float(a[2]),
            'weight_d':float(a[3]),
            'weight_phi':float(a[4]),
            'spread_threshold':float(a[5]),
            'burn_speed':float(a[6]),
            'spot_sigma':float(a[7]),
            'brand_sigma':float(a[8])
        })
        
        CA.run()
        
        np.savetxt(a[1], CA.CA_burnt_pixels, delimiter=',')
    
    else:
        CA.plotting = True
        CA.consts['gif_dir'] = a[1]
        CA.run()
    
if __name__=='__main__':
    main()
