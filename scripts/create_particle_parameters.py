'''
Created on Jan 28, 2021

@author: Stephen
'''
import numpy as np
from schedula.utils.cst import SELF

number_of_decimals = "11"
np.set_printoptions(formatter={'float': lambda x: "{0:0." + number_of_decimals + "f}".format(x)+"0"})

zero = [0.9999999999984769, 1.7453292519934434e-6, 0.0]

class Particles():
    def __init__(self):
        
        self.N = []
        self.masses = []
        self.Z = []
        self.E = []
        self.Ec = []
        self.Es = []
        self.positions = []
        self.directions = []
        self.interaction_index = []
        
    
    def add_particle_species(self, number_of_parts, N,  mass, Z, E, Ec, Es, position = None, direction = None, interaction_index = None):
        for _ in range(number_of_parts):
            self.N.append(N)
            self.masses.append(mass)
            self.Z.append(Z)
            self.E.append(E)
            self.Ec.append(Ec)
            self.Es.append(Es)
            
        
        
            if position != None:
                if type(position) == list:
                    self.positions.append(position)
            else:
                self.positions.append([0.0, 0.0, 0.0])
            if direction != None:
                if type(direction) == list:
                    self.directions.append(np.asarray(direction))
            else:
                self.directions.append(zero)
            
            if interaction_index != None:
                if type(interaction_index) == list:
                    self.interaction_index.append(interaction_index)
            else:
                self.interaction_index.append(0)
        
        
        return True
    
    def write_to_file(self, dump_to_file = False):
        
        import decimal
        decimal.getcontext().prec = int(number_of_decimals) + 2
        def have_ending_zeros(lis):
            #lis = np.asarray(lis)
            
            if type(lis[0]) == tuple:
                for i in range(len(lis)):
                    lis[i] = list(lis[i])
            
            for i in range(len(lis)):
                if type(lis[i]) == np.float64 or type(lis[i]) == np.float32 or type(lis[i]) == int or type(lis[i]) == float:
                    lis[i] =  decimal.Decimal(("{0:0." + str(len(str(lis[i]))) + "f}").format(lis[i])+"0")#.quantize(Decimal("1." + "0"*(int(number_of_decimals))))
                    #lis[i] =  decimal.Decimal(lis[i]).quantize(Decimal("1." + "0"*(int(number_of_decimals))))
                else:
                    if type(lis[i]) == type(decimal.Decimal(1.0)):
                                pass
                    for j in range(len(lis[i])):
                        if type(lis[i][j]) == np.float64 or type(lis[i][j]) == np.float32 or type(lis[i][j]) == int or type(lis[i][j]) == float or True:
                            lis[i][j] =  decimal.Decimal(("{0:0." + str(len(str(lis[i][j]))) + "f}").format(lis[i][j])+"0")#.quantize(decimal.Decimal("1." + "0"*(int(number_of_decimals))))
                            #lis[i][j] =  decimal.Decimal(lis[i][j]).quantize(decimal.Decimal("1." + "0"*(int(number_of_decimals))))
                        else:
                            if type(lis[i][j]) == type(decimal.Decimal(1.0)):
                                pass
                            for k in range(len(lis[i][j])):
                                if type(lis[i][j][k]) == np.float64 or type(lis[i][j][k]) == np.float32 or type(lis[i][j][k]) == int or type(lis[i][j][k]) == float:
                                    lis[i][j][k] =  decimal.Decimal(("{0:0." + str(len(str(lis[i][j][k]))) + "f}").format(lis[i][j][k])+"0")#.quantize(Decimal("1." + "0"*(int(number_of_decimals))))
                                    #lis[i][j][k] =  decimal.Decimal(lis[i][j][k]).quantize(Decimal("1." + "0"*(int(number_of_decimals))))
                                else:
                                    pass

            return lis
        
        file = open("Particles.toml", "w")
        file.seek(0)
        temp_dict = {
            "particle_parameters" :{
                "length_unit" : "MICRON",
                "energy_unit" : "EV",
                "mass_unit" : "AMU",
                "N" : [int(i) for i in self.N],
                "m" : have_ending_zeros(self.masses),
                "Z" : have_ending_zeros(self.Z),
                "E" : have_ending_zeros(self.E),
                "Ec" : have_ending_zeros(self.Ec),
                "Es" : have_ending_zeros(self.Es),
                "pos" : have_ending_zeros(self.positions),
                "dir" : have_ending_zeros(self.directions),
                "interaction_index" : [int(i) for i in self.interaction_index],
                "particle_input_filename" : ""
                }
            }
        if dump_to_file:
            import toml
            toml.dump(temp_dict, file)
        return temp_dict

if __name__ == "__main__":
    import random
    particle = Particles()
    
    for i in range(10):
        particle.add_particle_species(1, 1, 1.008, 1, 100000.0, .5, 10.0, position = [0.0,random.randint(0,50.0),0.0])
    particle.write_to_file(True)