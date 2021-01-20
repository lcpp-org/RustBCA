'''
Created on Jan 16, 2021

@author: Stephen
'''
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri
import scipy.spatial
import random
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


number_of_decimals = "10"
np.set_printoptions(formatter={'float': lambda x: "{0:0." + number_of_decimals + "f}".format(x)})


class Mesh():
    """
    This is the class that contains all the stuff to create a mesh for RustBCA. currently the only default shapes are N-point polygons and random points within the simulation boundaries
    """
    def __init__(self, length_unit, xmax, xmin, ymax, ymin, energy_barrier_thickness=1.7645653881793786e-4):
        """
        returns none.
        The init functions requires the length units, x limits and y limits of the simulation.
        The x and y limits are assumed to be on the points of a square.
        """
        
        self.trianglelist = []
        
        if not type(xmax) == int:
            print("xmax should be a int")
            assert type(xmax) == int
        if not type(xmin) == int:
            print("xmin should be a int")
            assert type(xmin) == int
        if not type(ymax) == int:
            print("ymax should be a int")
            assert type(ymax) == int
        if not type(ymin) == int:
            print("ymin should be a int")
            assert type(ymin) == int
        
        
        self.Simulation_boundaries = [[xmax, ymax], [xmax, ymin], [xmin, ymax], [xmin, ymin]]
        self.points = []
        self.shapes = []
        self.xmax, self.xmin, self.ymax, self.ymin = xmax, xmin, ymax,ymin
        if not type(length_unit) == str:
            print("length_units should be a string")
            assert type(length_unit) == str
        self.length_unit = length_unit
        self.electronic_stopping_corrections = []
        self.energy_barrier_thickness = energy_barrier_thickness
    
    def N_gon(self, radius, n_points, number_densities, x_offset = 0.0, y_offset = 0.0, theta_offset = 0.0):
        """
        returns true on completion
        This creates a polygon with n-points.
        The x and y offset determine the center point
        theta_offset is in radians
        
        easy shape examples:
        circle: N_gon(5,100,1)
        square: N_gon(2, 4, 2)
        offset square: N_gon(5, 4, 2, 1, 1, np.pi/4)  
        """
        
        n_points -= 0
        dtheta = np.pi*2/n_points
        
        temp_points = []
        for i in range(n_points):
            temp_points.append(
                Point(x_offset + radius*math.cos(dtheta*i + theta_offset), y_offset + radius*math.sin(dtheta*i + theta_offset))
                )

        #print(temp_points)
        poly = Polygon(temp_points)
        #print(poly)
        center = Point(x_offset,y_offset)
        #print(zero.within(poly))
        #print(poly.covers(zero))
        if center.within(poly):
            temp_points.append(center)
        #poly = Polygon(temp_points)
        #print(zero.within(poly))
        
        self.points += temp_points
        #temp_points = np.asarray(temp_points)
        
        self.shapes.append([poly, number_densities])
        
        return True
    
    def add_Uniform_random(self, n_points):
        """
        returns True on completion.
        Adds uniformly random points to break up the triangles into fairly random triangles.
        **Note**
        Will cause some shapes to have parts of them be labeled the incorrect number density. 
        """
        temp_points = []
        for _ in range(n_points):
            temp_points.append(
                Point(random.uniform(self.xmin, self.xmax), random.uniform(self.ymin, self.ymax))
                )
        #self.shapes.append(Polygon(temp_points))
        self.points += temp_points
        return True
    
    def get_Points(self):
        """
        Deprecated - Don't use
        returns all of the individual points in the simulation as a np array
        """
        temp = []
        for point in self.points:
            temp.append([point.x,point.y])
        
        return np.asarray(temp)
    
    def generate_Triangles(self):
        """
        returns a scipy.spatial.Delauny object
        Generates triangles via the Delauny method:
        see scipy.spatial.Delaunay
        """
        return scipy.spatial.Delaunay(self.get_Points())
    
    def return_Triangles(self):
        """
        returns 2 lists : points, material densities
        Creates and Correlates the triangles to their number densities.
        """
        points = self.get_Points()
        tri = self.generate_Triangles()
        
        point_output = []
        triangles = []
        
        for i in range(len(tri.simplices)):
            point_output.append([points[tri.simplices[i,0], 0], points[tri.simplices[i,1], 0], points[tri.simplices[i,2], 0], points[tri.simplices[i,0], 1], points[tri.simplices[i,1], 1], points[tri.simplices[i,2], 1]])
            triangles.append(Polygon([points[tri.simplices[i,0]], points[tri.simplices[i,1]], points[tri.simplices[i,2]]]))
            
        #print("Length of tri.simplices " + str(len(tri.simplices)))
        #print("Length of triangles " + str(len(triangles)))
        
        temp_material_densities = [0]*(len(triangles))
        self.electronic_stopping_corrections = [1.0]*(len(triangles))
        #print("Len of triangles " + str(len(triangles)))
        for i, triangle in enumerate(triangles):
            for shape, number_densities in self.shapes:
                if shape.contains(triangle):
                    temp_material_densities[i] = number_densities
        
        #print(point_output)
        return point_output, temp_material_densities
    
    def print_Triangles(self):
        """
        returns True on completion
        Use only as a check with number densities as integers from 0 to 3
        Plots the triangles with their colors based on their number densities
        Number densities should be integers currently
        """
        triangle_list, material_densities = self.return_Triangles()
    
        #print(len(triangle_list))
    
        #print(len(material_densities))
        
        self.color_dict = {0:"b", 1:"g", 2:"r", 3:"c", 4:"m"}
        
        triangles = []
        for t in triangle_list:
            triangles.append(
                matplotlib.tri.Triangulation(t[0:3], t[3:6])
                )
    
        #print(len(triangles))
        #print(type(triangles[0]))
        
        #plt.triplot(points[:,0], points[:,1], tri.simplices)
        plt.xlim(self.Simulation_boundaries[1][0], self.Simulation_boundaries[-1][0])
        plt.ylim(self.Simulation_boundaries[0:2][1])
        for i, triangle in enumerate(triangles):
            plt.triplot(triangle, self.color_dict[material_densities[i]]+ "-", )
        plt.show()
        return True
    
    def write_to_file(self, dump_to_file = False):
        '''
        returns a dictionary.
        Turn the required mesh2D stuff into a format that TOML and RustBCA like. 
        Optional to write to a .toml file instead of returning a dictionary
        '''
        triangle_list, material_densities = self.return_Triangles()
        Simulation_Boundaries = self.Simulation_boundaries
        
        material_boundary_points = []
        for shape, _ in self.shapes:
            material_boundary_points += shape.exterior.coords
        for i, point in reversed(list(enumerate(material_boundary_points))):
            for shape, _ in self.shapes:
                pointPoint = Point(point)
                if pointPoint.within(shape):
                    material_boundary_points.pop(i)
        
        import itertools
        material_boundary_points.sort()
        material_boundary_points = list(material_boundary_points for material_boundary_points,_ in itertools.groupby(material_boundary_points))
        #print(len(material_boundary_points))
        
        electronic_stopping_correction_factors = self.electronic_stopping_corrections
        
        file = open("Mesh2D.toml", "w")
        file.seek(0)
        
        temp_dict = {
            "[mesh_2d_input]" : {
                "length_unit":self.length_unit,
                "energy_barrier_thickness": self.energy_barrier_thickness,
                "triangles": triangle_list,
                "densities":material_densities,
                "material_boundary_points":material_boundary_points,
                "simulation_boundary_points":Simulation_Boundaries,
                "electronic_stopping_correction_factors":electronic_stopping_correction_factors
            }
        }
        
        if dump_to_file:
            import toml
            toml.dump(temp_dict, file)
        
        #print(type(triangle_list))
        
        return temp_dict
if __name__ == "__main__":
    
    mesh = Mesh("MICRON", 10,-10,10,-10)
    
    mesh.N_gon(5,100, [ 6.5E+10, 6.5E+10,])
    #mesh.N_gon(1.5, 4, 2, 5, 0, np.pi/4)
    #mesh.add_Uniform_random(10)
    #mesh.print_Triangles()
    mesh.write_to_file()
    
    
    #code to plot out the example triangles
    '''triangles = [ [ 0.0, 0.0, 0.0012558103905862675, 0.0, 0.02, 0.019960534568565433], [ 0.0, 0.0012558103905862675, 0.0025066646712860853, 0.0, 0.019960534568565433, 0.019842294026289557], [ 0.0, 0.0025066646712860853, 0.003747626291714493, 0.0, 0.019842294026289557, 0.019645745014573775], [ 0.0, 0.003747626291714493, 0.004973797743297096, 0.0, 0.019645745014573775, 0.01937166322257262], [ 0.0, 0.004973797743297096, 0.0061803398874989484, 0.0, 0.01937166322257262, 0.019021130325903073], [ 0.0, 0.0061803398874989484, 0.00736249105369356, 0.0, 0.019021130325903073, 0.018595529717765027], [ 0.0, 0.00736249105369356, 0.008515585831301454, 0.0, 0.018595529717765027, 0.01809654104932039], [ 0.0, 0.008515585831301454, 0.009635073482034308, 0.0, 0.01809654104932039, 0.01752613360087727], [ 0.0, 0.009635073482034308, 0.010716535899579934, 0.0, 0.01752613360087727, 0.0168865585100403], [ 0.0, 0.010716535899579934, 0.011755705045849463, 0.0, 0.0168865585100403, 0.016180339887498948], [ 0.0, 0.011755705045849463, 0.012748479794973795, 0.0, 0.016180339887498948, 0.015410264855515783], [ 0.0, 0.012748479794973795, 0.013690942118573775, 0.0, 0.015410264855515783, 0.014579372548428232], [ 0.0, 0.013690942118573775, 0.014579372548428232, 0.0, 0.014579372548428232, 0.013690942118573773], [ 0.0, 0.014579372548428232, 0.015410264855515785, 0.0, 0.013690942118573773, 0.012748479794973793], [ 0.0, 0.015410264855515785, 0.016180339887498948, 0.0, 0.012748479794973793, 0.011755705045849461], [ 0.0, 0.016180339887498948, 0.0168865585100403, 0.0, 0.011755705045849461, 0.01071653589957993], [ 0.0, 0.0168865585100403, 0.017526133600877274, 0.0, 0.01071653589957993, 0.009635073482034304], [ 0.0, 0.017526133600877274, 0.018096541049320392, 0.0, 0.009635073482034304, 0.008515585831301454], [ 0.0, 0.018096541049320392, 0.01859552971776503, 0.0, 0.008515585831301454, 0.007362491053693557], [ 0.0, 0.01859552971776503, 0.019021130325903073, 0.0, 0.007362491053693557, 0.006180339887498949], [ 0.0, 0.019021130325903073, 0.01937166322257262, 0.0, 0.006180339887498949, 0.004973797743297095], [ 0.0, 0.01937166322257262, 0.019645745014573775, 0.0, 0.004973797743297095, 0.0037476262917144902], [ 0.0, 0.019645745014573775, 0.019842294026289557, 0.0, 0.0037476262917144902, 0.0025066646712860853], [ 0.0, 0.019842294026289557, 0.019960534568565433, 0.0, 0.0025066646712860853, 0.001255810390586266], [ 0.0, 0.019960534568565433, 0.02, 0.0, 0.001255810390586266, -3.2162452993532734e-18], [ 0.0, 0.02, 0.019960534568565433, 0.0, -3.2162452993532734e-18, -0.001255810390586268], [ 0.0, 0.019960534568565433, 0.019842294026289557, 0.0, -0.001255810390586268, -0.0025066646712860875], [ 0.0, 0.019842294026289557, 0.019645745014573772, 0.0, -0.0025066646712860875, -0.0037476262917144967], [ 0.0, 0.019645745014573772, 0.01937166322257262, 0.0, -0.0037476262917144967, -0.004973797743297097], [ 0.0, 0.01937166322257262, 0.019021130325903073, 0.0, -0.004973797743297097, -0.006180339887498951], [ 0.0, 0.019021130325903073, 0.01859552971776503, 0.0, -0.006180339887498951, -0.00736249105369356], [ 0.0, 0.01859552971776503, 0.01809654104932039, 0.0, -0.00736249105369356, -0.008515585831301454], [ 0.0, 0.01809654104932039, 0.01752613360087727, 0.0, -0.008515585831301454, -0.00963507348203431], [ 0.0, 0.01752613360087727, 0.0168865585100403, 0.0, -0.00963507348203431, -0.010716535899579938], [ 0.0, 0.0168865585100403, 0.016180339887498948, 0.0, -0.010716535899579938, -0.011755705045849461], [ 0.0, 0.016180339887498948, 0.015410264855515785, 0.0, -0.011755705045849461, -0.012748479794973795], [ 0.0, 0.015410264855515785, 0.014579372548428228, 0.0, -0.012748479794973795, -0.013690942118573775], [ 0.0, 0.014579372548428228, 0.01369094211857377, 0.0, -0.013690942118573775, -0.014579372548428234], [ 0.0, 0.01369094211857377, 0.01274847979497379, 0.0, -0.014579372548428234, -0.015410264855515788], [ 0.0, 0.01274847979497379, 0.011755705045849465, 0.0, -0.015410264855515788, -0.016180339887498948], [ 0.0, 0.011755705045849465, 0.010716535899579934, 0.0, -0.016180339887498948, -0.0168865585100403], [ 0.0, 0.010716535899579934, 0.009635073482034304, 0.0, -0.0168865585100403, -0.01752613360087727], [ 0.0, 0.009635073482034304, 0.00851558583130145, 0.0, -0.01752613360087727, -0.018096541049320392], [ 0.0, 0.00851558583130145, 0.007362491053693555, 0.0, -0.018096541049320392, -0.01859552971776503], [ 0.0, 0.007362491053693555, 0.006180339887498942, 0.0, -0.01859552971776503, -0.019021130325903073], [ 0.0, 0.006180339887498942, 0.004973797743297097, 0.0, -0.019021130325903073, -0.01937166322257262], [ 0.0, 0.004973797743297097, 0.0037476262917144915, 0.0, -0.01937166322257262, -0.019645745014573775], [ 0.0, 0.0037476262917144915, 0.002506664671286082, 0.0, -0.019645745014573775, -0.019842294026289557], [ 0.0, 0.002506664671286082, 0.0012558103905862628, 0.0, -0.019842294026289557, -0.019960534568565433], [ 0.0, 0.0012558103905862628, -6.432490598706547e-18, 0.0, -0.019960534568565433, -0.02], [ 0.0, -6.432490598706547e-18, -0.0012558103905862669, 0.0, -0.02, -0.019960534568565433], [ 0.0, -0.0012558103905862669, -0.0025066646712860858, 0.0, -0.019960534568565433, -0.019842294026289557], [ 0.0, -0.0025066646712860858, -0.0037476262917144954, 0.0, -0.019842294026289557, -0.019645745014573772], [ 0.0, -0.0037476262917144954, -0.0049737977432971, 0.0, -0.019645745014573772, -0.01937166322257262], [ 0.0, -0.0049737977432971, -0.0061803398874989545, 0.0, -0.01937166322257262, -0.019021130325903073], [ 0.0, -0.0061803398874989545, -0.007362491053693567, 0.0, -0.019021130325903073, -0.018595529717765024], [ 0.0, -0.007362491053693567, -0.008515585831301454, 0.0, -0.018595529717765024, -0.01809654104932039], [ 0.0, -0.008515585831301454, -0.009635073482034308, 0.0, -0.01809654104932039, -0.01752613360087727], [ 0.0, -0.009635073482034308, -0.010716535899579936, 0.0, -0.01752613360087727, -0.0168865585100403], [ 0.0, -0.010716535899579936, -0.011755705045849468, 0.0, -0.0168865585100403, -0.016180339887498944], [ 0.0, -0.011755705045849468, -0.0127484797949738, 0.0, -0.016180339887498944, -0.015410264855515781], [ 0.0, -0.0127484797949738, -0.013690942118573775, 0.0, -0.015410264855515781, -0.014579372548428232], [ 0.0, -0.013690942118573775, -0.014579372548428232, 0.0, -0.014579372548428232, -0.013690942118573773], [ 0.0, -0.014579372548428232, -0.015410264855515788, 0.0, -0.013690942118573773, -0.01274847979497379], [ 0.0, -0.015410264855515788, -0.016180339887498948, 0.0, -0.01274847979497379, -0.011755705045849465], [ 0.0, -0.016180339887498948, -0.016886558510040308, 0.0, -0.011755705045849465, -0.010716535899579927], [ 0.0, -0.016886558510040308, -0.01752613360087727, 0.0, -0.010716535899579927, -0.009635073482034306], [ 0.0, -0.01752613360087727, -0.018096541049320396, 0.0, -0.009635073482034306, -0.008515585831301443], [ 0.0, -0.018096541049320396, -0.01859552971776503, 0.0, -0.008515585831301443, -0.007362491053693556], [ 0.0, -0.01859552971776503, -0.019021130325903073, 0.0, -0.007362491053693556, -0.006180339887498951], [ 0.0, -0.019021130325903073, -0.019371663222572624, 0.0, -0.006180339887498951, -0.004973797743297089], [ 0.0, -0.019371663222572624, -0.019645745014573775, 0.0, -0.004973797743297089, -0.003747626291714493], [ 0.0, -0.019645745014573775, -0.019842294026289557, 0.0, -0.003747626291714493, -0.0025066646712860745], [ 0.0, -0.019842294026289557, -0.019960534568565433, 0.0, -0.0025066646712860745, -0.001255810390586264], [ 0.0, -0.019960534568565433, -0.02, 0.0, -0.001255810390586264, -3.673940397442059e-18], [ 0.0, -0.02, -0.019960534568565433, 0.0, -3.673940397442059e-18, 0.0012558103905862745], [ 0.0, -0.019960534568565433, -0.019842294026289557, 0.0, 0.0012558103905862745, 0.0025066646712860845], [ 0.0, -0.019842294026289557, -0.019645745014573772, 0.0, 0.0025066646712860845, 0.003747626291714503], [ 0.0, -0.019645745014573772, -0.01937166322257262, 0.0, 0.003747626291714503, 0.004973797743297099], [ 0.0, -0.01937166322257262, -0.019021130325903073, 0.0, 0.004973797743297099, 0.006180339887498945], [ 0.0, -0.019021130325903073, -0.018595529717765024, 0.0, 0.006180339887498945, 0.007362491053693565], [ 0.0, -0.018595529717765024, -0.018096541049320392, 0.0, 0.007362491053693565, 0.008515585831301452], [ 0.0, -0.018096541049320392, -0.017526133600877267, 0.0, 0.008515585831301452, 0.009635073482034314], [ 0.0, -0.017526133600877267, -0.0168865585100403, 0.0, 0.009635073482034314, 0.010716535899579936], [ 0.0, -0.0168865585100403, -0.01618033988749894, 0.0, 0.010716535899579936, 0.011755705045849473], [ 0.0, -0.01618033988749894, -0.015410264855515781, 0.0, 0.011755705045849473, 0.0127484797949738], [ 0.0, -0.015410264855515781, -0.014579372548428232, 0.0, 0.0127484797949738, 0.013690942118573773], [ 0.0, -0.014579372548428232, -0.013690942118573766, 0.0, 0.013690942118573773, 0.014579372548428239], [ 0.0, -0.013690942118573766, -0.012748479794973793, 0.0, 0.014579372548428239, 0.015410264855515788], [ 0.0, -0.012748479794973793, -0.011755705045849453, 0.0, 0.015410264855515788, 0.016180339887498955], [ 0.0, -0.011755705045849453, -0.010716535899579927, 0.0, 0.016180339887498955, 0.016886558510040305], [ 0.0, -0.010716535899579927, -0.009635073482034308, 0.0, 0.016886558510040305, 0.01752613360087727], [ 0.0, -0.009635073482034308, -0.008515585831301445, 0.0, 0.01752613360087727, 0.018096541049320396], [ 0.0, -0.008515585831301445, -0.007362491053693557, 0.0, 0.018096541049320396, 0.01859552971776503], [ 0.0, -0.007362491053693557, -0.006180339887498935, 0.0, 0.01859552971776503, 0.019021130325903076], [ 0.0, -0.006180339887498935, -0.00497379774329709, 0.0, 0.019021130325903076, 0.019371663222572624], [ 0.0, -0.00497379774329709, -0.0037476262917144937, 0.0, 0.019371663222572624, 0.019645745014573775], [ 0.0, -0.0037476262917144937, -0.002506664671286076, 0.0, 0.019645745014573775, 0.019842294026289557], [ 0.0, -0.002506664671286076, -0.0012558103905862654, 0.0, 0.019842294026289557, 0.019960534568565433], [ 0.0, -0.0012558103905862654, -4.898587196589413e-18, 0.0, 0.019960534568565433, 0.02]]
    triangles = np.asarray(triangles)
    print(np.shape(triangles))
    
    points = []
    
    for i in range(len(triangles)):
        points.append([triangles[i,0], triangles[i,3]])
        points.append([triangles[i,1], triangles[i,4]])
        points.append([triangles[i,2], triangles[i,5]])
    
    
    import itertools
    points.sort()
    points = list(points for points,_ in itertools.groupby(points))
    
    points = np.asarray(points)
    print(np.shape(points))
    plt.triplot(points[:,0], points[:,1])
    plt.show()'''