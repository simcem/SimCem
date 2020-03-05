from simcem import *
import math

class Kiln:
    def __init__(self, innerRadius, layers,L):
        self.innerRadius = innerRadius
        self.insulationLayers = layers
        self.length = L
        
    def finish(self):
        CS_list = ["T1","T2","T3","T4"] #Coarse sand
        FS_list = ["T5","T6","T7","T8","T9"] #Fine Sand
        #EXPERIMENTS 10-13 are on Pet Coke, but missing!
        PC_list = ["T10","T11","T12","T13"] #Pet Coke
        #EXPERIMENTS 14-23 are on Limestone but missing!
        LS_list = ["T14","T15","T16","T17","T18","T19","T20","T21","T22","T23"] #Limestone

        self.k_expr = Expr("2.5")
        
        if self.kilnName == "Tscheng":
            #Not sure which Tscheng experiments used limestone, so all fitted to (Ottawa) sand
            solidComponents={'SiO2':1.0}
            self.bulk_density = 1650 #Tscheng thesis pdf pg.80
            self.solid_density = 2627 #Tscheng thesis pdf pg.80
        if self.kilnName == "KDO":
            solidComponents={'CaCO3':1.0}
            self.bulk_density = 2000
            self.solid_density = 2900
        if self.ID in CS_list:
            solidComponents={'SiO2':1.0}
            self.bulk_density = 1460  #Barr thesis pdf pg.98
            self.solid_density = 2627 #Barr thesis pdf pg.98
        if self.ID in FS_list:
            solidComponents={'SiO2':1.0}
            self.bulk_density = 1520  #Barr thesis pdf pg.98
            self.solid_density = 2627 #Barr thesis pdf pg.98
        if self.ID in LS_list:
            solidComponents={'CaCO3':1.0},
            self.bulk_density = 1680 #Barr thesis pdf pg.98, not correct when calcined!
            self.solid_density = 2300 #Assumed from internet searches!

        self.solidIn = Components(solidComponents) * (self.solidFeedRate / 1000.0)
        
class BarrKiln(Kiln):
    def __init__(self,solidFeedRate,primaryAir,secondaryAir,inGasFlow,solidLoading,RPM,particleDiameter):
        layers = [
            ["RefractoryBrick", 0.093, "0.2475 * (1.0 + 5.85e-4* T)"],#Refractory brick
            ["Steel", 0.006, "57"]#Shell
        ]
        Kiln.__init__(self, innerRadius=0.2055, layers=layers, L = 5.5)
        self.kilnName = "Barr"
        self.inGasFlow = inGasFlow#L/s
        self.solidFeedRate = solidFeedRate
        self.primaryAir = primaryAir#L/s
        self.secondaryAir = secondaryAir#L/s
        self.solidLoading = solidLoading #%
        self.RPM = RPM
        self.particleDiameter = particleDiameter

        
class BarrKilnT4(BarrKiln):
    def __init__(self):
        BarrKiln.__init__(self, solidFeedRate=62.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 17.4, secondaryAir = 43.0,inGasFlow = 1.97, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.0025#m
        )
        self.bedEmissivity = 0.9
        self.shellEmissivity = 0.80 
        self.wallEmissivity = 0.85   
        self.angleOfRepose = 31.0#deg
        self.incline = 0.0
        self.ID = "T4"
        self.y_gas_off_wall = [817.842,890.041,964.73,992.116,1007.05,1031.95,1056.85,1059.34,1081.74]
        self.x_gas_off_wall = [0.101266,0.898734,2.16456,2.51899,2.8481,3.22785,3.92405,4.4557,4.97468]
        self.y_bed = [486.722,646.058,740.664,802.905,830.29,855.187,870.124,927.386,962.241,994.606]
        self.x_bed =  [0.101266,0.886076,1.43038,2.13924,2.50633,2.86076,3.20253,3.93671,4.50633,4.97468]
        self.y_gas_off_bed= [765.56,865.145,944.813,974.689,992.116,1017.01,1049.38,1081.74,1129.05]
        self.x_gas_off_bed = [0.113924,0.911392,2.16456,2.53165,2.87342,3.24051,3.94937,4.48101,4.97468]
        self.y_wall = [730.705,812.863,837.759,857.676,875.104,947.303,984.647]
        self.x_wall = [1.35443,2.31646,2.68354,3.03797,3.39241,4.41772,4.73418]
        #self.phase = "1a-qz"     
    
    
class BarrKilnT8(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=65.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.8, secondaryAir = 40.01,inGasFlow = 2.00, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.00058#m
                          )	
        self.ID = "T8"
        self.y_gas_off_wall = [834.648,876.944,959.68,974.382,1010.67,1031.86,1059.12,1105.45]
        self.x_gas_off_wall =  [0.113848,0.884104,2.14718,2.51063,2.8539,3.19637,3.91249,4.94962]
        self.y_bed = [450.82,618.199,712.463,780.714,810.51,823.068,857.199,903.877,950.693,995.478]
        self.x_bed =  [0.11471,0.870304,1.44066,2.13764,2.50189,2.85456,3.19772,3.90421,4.4827,4.94376]
        self.y_gas_off_bed= [741.92,838.121,935.95,963.601,989.108,1046.18,1120.55]
        self.x_gas_off_bed = [0.119569,0.8927,2.15658,2.51006,2.85275,3.92246,4.93976]
        self.y_wall = [710.422,791.3,818.951,837.977,865.629,912.582,987.026]
        self.x_wall = [1.33389,2.31954,2.67301,3.02603,3.3795,3.83,4.78331]  
        
        
        
class BarrKilnT2(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=62.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 16.5, secondaryAir = 40.6,inGasFlow = 1.02, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.0025#m
                          )	
        self.ID = "T2"
        self.y_gas_off_wall = [612.5,660.714,707.143,726.786,726.786,744.643,757.143,771.429,787.5]
        self.x_gas_off_wall = [0.111173,0.883217,2.15076,2.50561,2.85074,3.20552,3.89633,4.41024,4.93308]
        self.y_bed = [371.429,528.571,591.071,650.0,673.214,700.0,730.357]
        self.x_bed = [0.109355,1.43489,2.4996,3.90929,4.46784,4.93805,5.490]
        self.y_gas_off_bed= [589.286,619.643,687.5,707.143,716.071,757.143,764.286,783.929]
        self.x_gas_off_bed = [0.110145,0.89025,2.13219,2.50474,2.85027,3.90518,4.42762,4.94177]
        self.y_wall = [535.714,583.929,601.786,616.071,628.571,651.786,687.5,712.5]
        self.x_wall = [1.32901,2.33115,2.66822,3.01399,3.38622,3.82088,4.38883,4.77046]         
        
class BarrKilnT3(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=62.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 17.4, secondaryAir = 40.6,inGasFlow = 1.42, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.0025#m
                          )	
        self.ID = "T3"
        self.y_gas_off_wall = [701.567,765.204,824.451,837.618,848.589,855.172,890.282,907.837]
        self.x_gas_off_wall = [0.120382,0.899883,2.15441,2.51119,2.8571,3.19208,3.9059,4.93236]
        self.y_bed = [464.577,574.295,642.32,688.401,712.539,725.705,749.843,778.37,806.897,857.367,815.674,852.978]
        self.x_bed = [0.10227,0.872394,1.43605,2.13941,2.49653,2.85331,3.19962,3.91325,4.48648,4.9416,5.25349,5.49222]
        self.y_gas_off_bed= [673.041,730.094,795.925,815.674,833.229,848.589,894.671]
        self.x_gas_off_bed = [0.119501,0.888,2.14273,2.51051,2.84583,3.20267,4.44599]
        self.y_wall = [626.959,694.984,714.734,732.288,749.843,789.342,822.257,850.784]
        self.x_wall =  [1.32758,2.3124,2.66938,3.0155,3.37241,3.83799,4.38976,4.77941]

class BarrKilnT5(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=58.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 9.4, secondaryAir = 19.8,inGasFlow = 0.68, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.00058 #m
                          )	
        self.ID = "T5"
        self.y_gas_off_wall = [585.573,640.609,714.213,735.246,749.782,770.839,808.566,818.314,836.834]
        self.x_gas_off_wall = [0.112579,0.887618,2.15241,2.51222,2.84976,3.19874,3.90714,4.43947,4.95091]
        self.y_bed = [400.107,453.009,509.162,538.921,553.41,587.555,640.552,685.131,701.588,735.828,768.005]
        self.x_bed = [0.117508,0.870662,2.13387,2.49448,2.85371,3.20386,3.91364,4.48166,4.93868,5.24547,5.49783]
        self.y_gas_off_bed= [557.191,570.781,613.845,643.628,675.591,727.188,778.003,800.863,808.476]
        self.x_gas_off_bed = [0.12086,0.892153,2.15418,2.50394,2.8539,3.20564,3.91522,4.43789,4.94834]
        self.y_wall = [504.345,563.319,584.375,603.227,624.284,662.579,694.069,732.483]
        self.x_wall = [1.34168,2.3123,2.66128,3.0209,3.36987,3.81802,4.38486,4.77879]

class BarrKilnT6(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=62.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 14.2, secondaryAir = 29.3,inGasFlow = 0.90, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.00058 #m
                          )	
        self.ID = "T6"
        self.y_gas_off_wall = [673.041,690.596,747.649,756.426,776.176,780.564,813.48,802.508,802.508]
        self.x_gas_off_wall = [0.114705,0.89474,2.16362,2.49953,2.85758,3.20411,3.91992,4.44971,4.96919]
        self.y_bed = [365.831,484.326,541.379,585.266,605.016,618.182,640.125,668.652,694.984,727.9]
        self.x_bed =  [0.111278,0.87435,1.43977,2.14527,2.50332,2.86108,3.20842,3.9132,4.48802,4.95491]
        self.y_gas_off_bed= [585.266,642.32,727.9,747.649,760.815,795.925,787.147,787.147]
        self.x_gas_off_bed = [0.121456,0.892501,2.1627,2.50994,2.86769,3.91911,4.449,4.95766]
        self.y_wall = [565.517,620.376,637.931,653.292,662.069,694.984,714.734,738.871]
        self.x_wall = [1.33266,2.32005,2.67801,3.03587,3.39342,3.83867,4.40235,4.79308]

class BarrKilnT7(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=63.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.4, secondaryAir = 43.0,inGasFlow = 1.04, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.00058 #m
                          )	
        self.ID = "T7"
        self.y_gas_off_wall = [651.066,659.607,703.487,709.999,716.512,725.246,742.71,742.469,755.581]
        self.x_gas_off_wall = [0.101999,0.88328,2.15416,2.50152,2.84887,3.19632,3.90207,4.44436,4.94384]
        self.y_bed = [428.843,504.06,530.475,550.161,572.224,587.62,605.248,629.374,660.229,662.244,684.235] 
        self.x_bed =   [0.103204,0.865686,1.43082,2.13667,2.49554,2.85413,3.19112,3.908,4.48417,4.95064,5.4722]
        self.y_gas_off_bed= [582.173,628.497,690.144,701.106,705.43,729.372,744.47]
        self.x_gas_off_bed =  [0.109856,0.88193,2.17528,2.51198,2.78332,3.91234,4.94336] 
        self.y_wall = [534.963,581.195,598.813,609.765,618.5,644.968,684.538]
        self.x_wall = [1.3334,2.31155,2.67023,3.02862,3.37607,3.82191,4.78892]

class BarrKilnT1(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=62.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 9.4, secondaryAir = 18.8,inGasFlow = 0.83, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.0025 #m
                          )	
        self.ID = "T1"
        self.y_gas_off_wall = [592.996,661.964,749.991,769.207,803.748,818.617,872.36,943.188]
        self.x_gas_off_wall = [0.0986207,0.884875,2.14983,2.50924,2.85898,3.19636,3.91636,4.95221]
        self.y_bed =  [400.487,454.172,514.676,555.323,606.878,751.305,809.75,854.971]
        self.x_bed =   [0.105399,0.868777,1.43746,2.13475,2.85458,4.48182,4.95272,5.49853]
        self.y_gas_off_bed= [564.547,576.646,629.676,662.016,692.196,744.25,848.3,891.362,927.877]
        self.x_gas_off_bed =  [0.107263,0.889111,2.15135,2.51178,2.85034,3.1906,3.91449,4.43844,4.95103] 
        self.y_wall = [530.129,587.89,609.293,635.099,652.128,704.025,764.543,812.15]
        self.x_wall = [1.33019,2.32165,2.68123,3.01945,3.37869,3.83824,4.39608,4.79022]
        
class BarrKilnT9(BarrKilnT4):
    def __init__(self):
        BarrKilnT4.__init__(self)
        BarrKiln.__init__(self, solidFeedRate=64.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.8, secondaryAir = 43.0,inGasFlow = 2.53, #L/s
                          solidLoading = 0.12,RPM = 1.5,particleDiameter = 0.00058 #m
                          )	
        self.ID = "T9"
        self.y_gas_off_wall = [874.748,966.177,1072.98,1090.49,1116.66,1136.34,1173.51,1219.5]
        self.x_gas_off_wall = [0.114226,0.908667,2.15839,2.51638,2.86396,3.22205,3.91647,4.96847]
        self.y_bed =  [503.907,701.579,797.252,882.137,921.324,951.841,978.025,1023.87,1076.17,1126.25]
        self.x_bed =  [0.118033,0.863472,1.44167,2.14921,2.49742,2.85603,3.22525,3.93092,4.49621,4.97481]
        self.y_gas_off_bed= [805.349,924.972,1044.8,1070.97,1173.51,1178.08,1252.04]
        self.x_gas_off_bed =  [0.110887,0.906685,2.17868,2.51544,3.9273,4.45782,4.97004]
        self.y_wall = [803.711,901.731,938.754,962.765,995.45,1041.19,1071.79,1113.17]
        self.x_wall =   [1.33375,2.32332,2.68224,3.04054,3.39925,3.84517,4.40942,4.80102]  

class TschengKiln(Kiln):
    def __init__(self,solidFeedRate,primaryAir,secondaryAir,solidLoading,RPM,particleDiameter):

        ##All page numbers refer to Tscheng's thesis

        #The thermocouples in the wall are placed between the
        #refactory cement and the steel shell (so are not equivalent
        #to our model's temperatures)

        #From 1.02-1.22m there is no insulation (due to the commutator
        #ring).  They recommend taking HT data from 1.25-1.78m only.
        #
        #There's also six 6.7cm x 10cm rollers (either 40cm or 60cm of
        #kiln length in total) which must be uninsulated, and end
        #boxes which have ~25cm of kiln within them (in total) which
        #are probably uninsulated.
        #
        #There's also the friction belt zone which is probably
        #uninsulated.
        #
        #In total, at least 85cm-95cm of uninsulated kiln, and the
        #whole kiln is 2.44m (>35%)!

        #There is also some confusion about the actual insulation
        #itself.  On pdf pg.87 of Tscheng's thesis, they say that
        #3.2mm asbestos cloth and 25.4mm fibreglass was used over the
        #ceramic paper; however, later the fibreglass was replaced
        #with 51mm of fibred asbestos pipe insulation due to
        #significant heat loss, reducing it from 50% to 20% of the
        #total gas enthalpic change. Both of these descriptions are
        #different to what is described earlier in the thesis and what
        #is in the paper (which are in agreement)

        #In all, this makes a difficult task to actually calculate the
        #heat lost. Here, we follow the sample calculations (pg 253
        #onward) to figure out that q_sw=57.2 W/m and q_gw = 97.2 W/m
        #leading to a total flux of q_loss=154.4 W/m. They are looking
        #at test A16, Z=1.25m, and interpolate the wall temperature as
        #T_w=369.2 K, solid temperature T_s=374.0 K. Thus the
        #effective UA is 2.17 W/m/K. As there is a nonlinear
        #relationship in our model (thanks to radiation on the outer
        #surface), we can only choose to match either UA or
        #q_loss. Taking q_loss, an effective thickness of fiberglass
        #is 14.87 mm as used below
        layers = [
            ["RefractoryCement", 0.001, "0.294"], #Taken from pg255 or pg266
            ["Steel", 0.00635, "45.2"],#Shell WARNING value thermal conductivity assumed
            ["Ceramic paper", 0.00635, "0.16"],#Ceramic paper (googled
            #["Fiberglass", 0.051, "0.08"],# Fiber glass (googled)
            ["Fiberglass", 0.01487, "0.08"],#
        ]
        Kiln.__init__(self, innerRadius=0.09425,#0.09525 according to TScheng on pdf pg.266.
                      layers=layers,
                      L = 2.44)
        self.kilnName = "Tscheng"
        self.solidFeedRate = solidFeedRate
        self.primaryAir = primaryAir / 1.225 / 3.6#convert from kg/hr to L/s
        self.secondaryAir = secondaryAir/ 1.225 / 3.6#convert from kg/hr to L/s
        self.inGasFlow = 0#L/s //This kiln only uses air!
        self.solidLoading = solidLoading #%
        self.RPM = RPM
        self.omega =  self.RPM*2*math.pi/60.0 #Kiln rotational velocity in Radians per seconds
        self.particleDiameter = particleDiameter
        
        #Create solid stream

class TschengKilnA11(TschengKiln):
    def __init__(self):
        TschengKiln.__init__(self, solidFeedRate=25.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6 , secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 3.0,particleDiameter =0.00073#m
                          )
        self.bedEmissivity = 0.9 
        self.shellEmissivity = 0.80 
        self.wallEmissivity = 0.85   
        self.angleOfRepose = 27.0#deg
        self.incline = 1.2
        self.ID = "A11"
        self.y_gas_off_wall = [450.0,486.0,524.0,574.0,635.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed = [334.0,356.0,378.0,431.0,521.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [320.0,351.0,397.0,500.0]
        self.x_wall = [0.31,0.91,1.52,2.13]
        #self.phase = "1a-qz"     
        
class TschengKilnA12(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate=25.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 3.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A12"
        self.incline = 1.2
        self.y_gas_off_wall = [425.0,460,489,535,592]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [323,339,356,397,473]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [314,338,372,445]
        self.x_wall =   [0.31,0.91,1.52,2.13]  
        
class TschengKilnA13(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate=25.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 3.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A13"
        self.incline = 1.2
        self.y_gas_off_wall = [402,423,457,493,538]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [321,337,354,392,447]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [313,333,368,429]
        self.x_wall =   [0.31,0.91,1.52,2.13]        
        
        
class TschengKilnA14(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 25.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 3.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A14"
        self.incline = 1.2
        self.y_gas_off_wall = [372.0,392,411,436,457]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [314,327,341,368,405]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [308,326,352,395]
        self.x_wall =   [0.31,0.91,1.52,2.13]         
        
class TschengKilnA15(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 14.2 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.50,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A15"
        self.incline = 1.2
        self.y_gas_off_wall = [330.6,340.0,348.3,358.8,372.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [306.7,308.9,312.8,326.7,348.3]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [301.0,306.0,316.0,333.0]
        self.x_wall =   [0.31,0.91,1.52,2.13]         
        
class TschengKilnA16(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 14.2 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.5,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A16"
        self.incline = 1.2
        self.y_gas_off_wall = [414.0,438.0,462.0,494.0,535.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [341.0,356.0,374.0,417.0,473.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [324.0,347.0,387.5,447.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
        
class TschengKilnA17(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 14.2 * 1000.0 / 3600.0,#g/s
                          primaryAir= 24.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.5,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A17"
        self.incline = 1.2
        self.y_gas_off_wall = [470.0,505.0,543.0,594.0,652.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [364.0,383.0,410.0,473.0,554.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [341.0,370.0,425.0,520.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

        
class TschengKilnA18(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 21.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.15,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A18"
        self.incline = 1.2
        self.y_gas_off_wall = [351.5,378.9,412.2,445.0,488.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [313.9,325.6,338.9,363.9,402.8]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [308.0,323.0,349.0,385.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

        
class TschengKilnA19(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 29.1 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.6,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A19"
        self.incline = 2.2
        self.y_gas_off_wall = [350.0,375.1,405.0,448.0,497.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [312.2,322.2,332.2,355.5,391.7]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [304.0,317.0,341.0,377.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

        
class TschengKilnA20(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 15.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.6,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A20"
        self.incline = 1.2
        self.y_gas_off_wall = [365.0,392.8,427.2,464.0,507.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [330.5,340.0,352.8,383.3,424.4]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [315.5,331.0,363.0,407.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA21(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 15.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.5,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A21"
        self.incline = 1.2
        self.y_gas_off_wall = [380.0,398.0,418.0,436.0,455.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [335.0,346.1,360.3,384.4,411.1]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [319.4,336.0,365.8,401.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
        
        
class TschengKilnA22(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 34.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A22"
        self.incline = 1.2
        self.y_gas_off_wall = [368.0,388.0,407.0,428.0,455.8]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [320.6,331.1,341.7,363.3,390.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [311.0,326.6,349.6,378.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
         
class TschengKilnA23(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 15.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.5,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A23"
        self.incline = 1.2
        self.y_gas_off_wall = [383.0,400.8,418.0,437.0,456.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [337.8,348.9,365.0,385.6,412.8]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [321.0,338.0,365.0,400.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA24(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 50.5 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 6.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A24"
        self.incline = 1.2
        self.y_gas_off_wall = [357.0,374.0,390.6,407.0,430.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [308.3,315.0,322.2,336.7,360.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [302.0,311.6,330.0,355.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA25(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 39.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.17,RPM = 1.50,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A25"
        self.incline = 3.40
        self.y_gas_off_wall = [361.0,377.2,391.5,410.0,430.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [305.6,313.6,322.2,338.9,362.1]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [300.0,313.4,331.0,360.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA26(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 34.6 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.2,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A26"
        self.incline = 2.2
        self.y_gas_off_wall = [370.0,385.0,398.9,414.0,437.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [313.8,323.3,333.3,350.0,376.1]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [306.0,321.6,340.0,364.8]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA27(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 34.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 50.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.2,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A27"
        self.incline = 2.2
        self.y_gas_off_wall = [369.0,380.0,388.9,401.0,413.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [316.7,326.1,334.4,351.7,373.9]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [306.8,323.0,341.6,366.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA28(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 52.7 * 1000.0 / 3600.0,#g/s
                          primaryAir= 50.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.10,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A28"
        self.incline = 3.00
        self.y_gas_off_wall = [361.0,371.0,381.1,395.0,410.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [312.8,318.9,327.9,341.7,358.9]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [304.0,315.5,333.0,353.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA29(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 19.4 * 1000.0 / 3600.0,#g/s
                          primaryAir= 50.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.10,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A29"
        self.incline = 1.2
        self.y_gas_off_wall = [376.0,386.0,396.0,406.0,417.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [328.9,339.4,353.9,370.6,391.1]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [320.0,338.8,350.0,380.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA30(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 18.2 * 1000.0 / 3600.0,#g/s
                          primaryAir= 50.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 1.60,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A30"
        self.incline = 2.2
        self.y_gas_off_wall = [378.0,388.0,399.0,410.0,420.5]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [329.4,341.1,357.2,373.9,394.4]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [320.0,339.0,360.0,387.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA31(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 66.3 * 1000.0 / 3600.0,#g/s
                          primaryAir= 50.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 6.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A31"
        self.incline = 2.2
        self.y_gas_off_wall = [353.0,365.0,375.6,387.5,404.5]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [308.3,313.9,321.5,331.7,350.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [302.3,313.6,326.0,345.6]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

class TschengKilnA32(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 81.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A32"
        self.incline = 2.00
        self.y_gas_off_wall = [407.0,417.8,425.1,437.5,450.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [333.3,350.6,366.7,391.1,417.3]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [321.1,347.4,375.5,407.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
                                                      
class TschengKilnA33(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 65.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A33"
        self.incline = 2.00
        self.y_gas_off_wall = [393.0,407.0,418.3,431.0,445.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [326.1,341.7,358.3,380.6,407.2]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [313.0,338.4,367.4,395.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA34(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 73.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A34"
        self.incline = 2.00
        self.y_gas_off_wall = [396.2,412.2,423.3,436.0,448.9]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [331.7,347.2,362.7,385.8,412.2]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [315.8,344.0,370.0,402.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA35(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 81.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A35"
        self.incline = 2.00
        self.y_gas_off_wall = [396.8,411.1,422.2,433.0,444.4]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [326.7,348.9,366.7,388.9,414.2]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [315.8,345.0,375.5,404.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

class TschengKilnA36(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A36"
        self.incline = 2.00
        self.y_gas_off_wall = [366.7,380.6,395.6,415.0,441.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [311.1,322.2,334.4,352.2,378.1]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [304.5,321.7,339.0,368.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
                                                
class TschengKilnA37(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A37"
        self.incline = 2.00
        self.y_gas_off_wall = [420.0,445.4,476.7,513.0,560.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [322.8,343.3,369.4,399.7,443.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [312.6,341.1,380.4,425.5]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
                                                     
class TschengKilnA38(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A38"
        self.incline = 2.00
        self.y_gas_off_wall = [395.0,417.0,440.6,471.0,505.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [318.9,334.4,355.0,380.0,415.3]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [312.4,335.4,361.0,401.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

class TschengKilnA39(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A39"
        self.incline = 2.00
        self.y_gas_off_wall = [420.6,446.5,475.0,512.0,559.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [323.9,343.3,369.4,400.0,445.4]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [312.0,344.0,380.0,425.5]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA40(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A40"
        self.incline = 2.00
        self.y_gas_off_wall = [385.6,417.0,461.1,510.0,560.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [318.3,331.7,352.0,377.8,428.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [311.0,333.0,364.0,413.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA41(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A41"
        self.incline = 2.00
        self.y_gas_off_wall = [358.3,376.0,410.6,447.0,505.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [308.9,317.8,332.2,348.3,387.2]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [300.4,317.6,335.0,372.8]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA42(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A42"
        self.incline = 2.00
        self.y_gas_off_wall = [351.7,360.0,384.4,413.5,460.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [306.1,312.2,321.7,336.1,366.7]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [299.0,312.6,327.0,356.8]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
        
class TschengKilnA43(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 18.6, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A43"
        self.incline = 2.00
        self.y_gas_off_wall = [375.6,400.0,440.0,478.0,535.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [313.3,325.0,345.0,365.0,412.2]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [306.0,325.0,352.0,397.5]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA44(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 36.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 50.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.11,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A44"
        self.incline = 2.00
        self.y_gas_off_wall = [398.0,416.7,433.9,452.5,475.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [323.9,338.9,360.0,384.2,419.4]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [312.0,338.0,371.0,407.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA45(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 12.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 65.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 0.9,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A45"
        self.incline = 2.00
        self.y_gas_off_wall = [412.0,425.0,434.4,441.1,448.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [361.1,373.9,391.7,416.7,434.4]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [339.8,365.0,400.0,428.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA46(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 13.3 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 1.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A46"
        self.incline = 2.00
        self.y_gas_off_wall = [382.1,397.2,413.3,427.3,444.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [333.3,345.6,361.7,386.7,411.1]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [314.0,338.4,369.4,400.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA47(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 35.8 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A47"
        self.incline = 2.00
        self.y_gas_off_wall = [367.8,382.2,398.9,416.0,441.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [312.8,322.8,336.1,352.8,380.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [306.0,323.6,344.2,370.5]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA48(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 35.8 * 1000.0 / 3600.0,#g/s
                          primaryAir= 65.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A48"
        self.incline = 2.00
        self.y_gas_off_wall = [400.0,408.0,420.6,432.0,447.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [327.8,341.1,360.0,382.0,409.4]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [318.0,340.2,369.0,395.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA49(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 11.7 * 1000.0 / 3600.0,#g/s
                          primaryAir= 65.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 0.9,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A49"
        self.incline = 2.00
        self.y_gas_off_wall = [419.0,427.0,435.0,442.3,448.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [368.9,378.3,397.8,419.4,438.3]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [349.6,374.4,402.5,430.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA50(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 35.8 * 1000.0 / 3600.0,#g/s
                          primaryAir= 95.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 3.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A50"
        self.incline = 2.00
        self.y_gas_off_wall = [404.0,410.2,417.2,424.0,431.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [375.0,377.8,391.7,407.2,418.9]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [355.0,374.0,395.0,415.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA51(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 15.8 * 1000.0 / 3600.0,#g/s
                          primaryAir= 95.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 1.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A51"
        self.incline = 2.00
        self.y_gas_off_wall = [406.6,412.2,417.6,422.0,427.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [370.6,377.2,390.6,402.7,418.3]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [351.4,372.2,393.3,413.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA52(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 11.3 * 1000.0 / 3600.0,#g/s
                          primaryAir= 34.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 1.0,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A52"
        self.incline = 2.00
        self.y_gas_off_wall = [384.0,398.1,416.7,428.0,443.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [335.5,348.9,368.3,391.7,418.9]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [319.0,341.5,373.0,405.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA53(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 16.1 * 1000.0 / 3600.0,#g/s
                          primaryAir= 95.5, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 1.00,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A53"
        self.incline = 2.00
        self.y_gas_off_wall = [403.3,409.5,414.4,413.5,422.8]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [374.4,375.6,388.9,404.4,415.6]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [347.0,369.0,395.0,410.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 
        
class TschengKilnA54(TschengKilnA11):
    def __init__(self):
        TschengKilnA11.__init__(self)
        TschengKiln.__init__(self, solidFeedRate= 12.0 * 1000.0 / 3600.0,#g/s
                          primaryAir= 81.0, secondaryAir = 0.0, #kg/hr
                          solidLoading = 0.065,RPM = 0.95,particleDiameter = 0.00073 #m
                          )	
        self.ID = "A54"
        self.incline = 2.00
        self.y_gas_off_wall = [419.0,427.0,434.4,442.0,450.0]
        self.x_gas_off_wall = [0.21,0.72,1.25,1.78,2.32]
        self.y_bed =  [373.3,381.7,397.2,418.9,435.0]
        self.x_bed =  [0.21,0.72,1.25,1.78,2.32]
        self.y_gas_off_bed=[]
        self.x_gas_off_bed=[]
        self.y_wall = [351.0,364.4,400.0,430.0]
        self.x_wall =   [0.31,0.91,1.52,2.13] 

class KDOKiln(Kiln):
    def __init__(self,solidFeedRate,primaryAir,secondaryAir,inGasFlow,solidLoading,RPM,particleDiameter,Tair = 298.15,Pair = 1.01325e5,Tgas = 298.15,Pgas = 1.01325e5):
      
        layers = [
            ["RefractoryBrick", 0.05, "1.39666667 + 2.5e-04 * T"],#Refractory brick
            ["SteelShell", 0.005, "57"]#Shell
            ]
        Kiln.__init__(self, innerRadius=0.15, layers=layers, L = 7.4)
        self.kilnName = "KDO"
        self.inGasFlow = inGasFlow#L/s
        self.solidFeedRate = solidFeedRate
        self.ambientTemperature = 303.15
        self.burner_radius = 0.05
        self.primaryAir = primaryAir#L/s
        self.secondaryAir = secondaryAir#L/s
        self.solidLoading = solidLoading #%
        self.RPM = RPM
        self.omega =  self.RPM*2*math.pi/60.0 #Kiln rotational velocity in Radians per seconds
        self.particleDiameter = particleDiameter
        
        #Use the adiabatic gas temperature for the outlet, assume room temperature solid inlet, and solve
        self.SolveMode = True

        
class KDOKilnTP21(KDOKiln):
    def __init__(self):
        KDOKiln.__init__(self, solidFeedRate=27.0 * 1000.0 / 3600.0,#g/s # Feed solid rate
                         primaryAir = 263 * 1000.0 / 3600.0, secondaryAir = 40 * 1000 / 3600 ,inGasFlow = 19 * 1000.0 / 3600.0, #L/s
                         solidLoading = 0.04,RPM = 1.8,particleDiameter = 0.0025#m
        )
        self.inO2Flow = 4 * 1000.0 / 3600 #L/s
        self.inSO2Flow = 23 / 60 #L/s
        self.bedEmissivity = 0.9 
        self.shellEmissivity = 0.80 
        self.wallEmissivity = 0.85   
        self.angleOfRepose = 31.0#deg
        self.incline = 0.0
        self.ID = "TP21"
        self.y_gas_off_wall = []#[1637.7] #I think these are adiabatic values from theo?
        self.x_gas_off_wall = []#[7.0]
        #self.y_bed = [532,730,829,1238,1250]
        #self.x_bed =  [0.19,1.69,3.32,6.06,6.51]
        #self.y_gas_off_bed= [1259]
        #self.x_gas_off_bed = [7.0]
        
        self.y_bed = [1259+273.15]
        self.x_bed =  [7.0]
        self.y_gas_off_bed= [642+273.15,532+273.15,730+273.15,829+273.15,1238+273.15,1250+273.15]
        self.x_gas_off_bed = [0.0,0.19,1.69,3.32,6.06,6.51]
        
        self.y_wall = [290+273.15,355+273.15,408+273.15,472+273.15,488+273.15,497+273.15]
        self.x_wall = [0.19,1.69,3.32,5.0,6.06,6.51]

class KDOKilnTP21r(KDOKiln):
    def __init__(self):
        KDOKiln.__init__(self, solidFeedRate=27.0 * 1000.0 / 3600.0,#g/s
                         primaryAir = 263 * 1000.0 / 3600.0, secondaryAir = 40 * 1000 / 3600, inGasFlow = 19 * 1000.0 / 3600.0, #L/s
                         solidLoading = 0.04,RPM = 1.8,particleDiameter = 0.0025#m
        )
        self.inO2Flow = 4 * 1000.0 / 3600 #L/s 
        self.inSO2Flow = 23 / 60 #L/s
        self.bedEmissivity = 0.9
        self.shellEmissivity = 0.80
        self.wallEmissivity = 0.85
        self.angleOfRepose = 31.0#deg
        self.incline = 0.0
        self.ID = "TP21r"
        self.y_gas_off_wall = []#[1440]
        self.x_gas_off_wall = []#[7.2]
        #self.y_bed = [532,730,829,1238,1250]
        #self.x_bed =  [0.19,1.69,3.32,6.06,6.51]
        #self.y_gas_off_bed= [1259]
        #self.x_gas_off_bed = [7.0]
        
        self.y_bed = [1259+273.15]
        self.x_bed =  [7.0]
        self.y_gas_off_bed= [532+273.15,730+273.15,829+273.15,1238+273.15,1250+273.15]
        self.x_gas_off_bed = [0.19,1.69,3.32,6.06,6.51]
        
        self.y_wall = [290+273.15,355+273.15,408+273.15,472+273.15,488+273.15,497+273.15]
        self.x_wall = [0.19,1.69,3.32,5.0,6.06,6.51]

class KDOKilnTP26(KDOKiln):
    def __init__(self):
        KDOKiln.__init__(self, solidFeedRate=9.65 * 1000.0 / 3600.0,#g/s
                          primaryAir= 70.0, secondaryAir = 0.0,inGasFlow = 4.9, #L/s
                          solidLoading = 0.08,RPM = 1.8,particleDiameter = 0.0025#m
                          )
        self.inO2Flow = 4 * 1000.0 / 3600 #L/s ##### WRONG!!!
        self.inSO2Flow = 23 / 60 #L/s ##### WRONG!!!
        self.bedEmissivity = 0.9 
        self.shellEmissivity = 0.80 
        self.wallEmissivity = 0.85   
        self.angleOfRepose = 31.0#deg
        self.incline = 0.0
        self.ID = "TP26"
        self.y_gas_off_wall = [1765.17]
        self.x_gas_off_wall =  [7.0]
        self.y_bed = [1280+273.15]
        self.x_bed =  [7.0]
        
        #self.y_gas_off_bed= [640,751]
        #self.x_gas_off_bed = [1.69,3.32]
        self.y_gas_off_bed= [597+273.15,503+273.15,640+273.15,751+273.15,1262+273.15]
        self.x_gas_off_bed = [0.0,0.19,1.69,3.32,6.51]
        
        self.y_wall = [275+273.15,340+273.15,352+273.15,381+273.15,386+273.15,391+273.15,414+273.15,415+273.15,463+273.15,440+273.15,498+273.15,490+273.15,443+273.15]
        self.x_wall = [0.19,1.69,1.97,2.53,2.83,3.32,3.7,4.09,5.0,5.3,6.06,6.51,7.08]
        #self.y_wall = [340,352,381,386,391,414,415]
        #self.x_wall = [1.69,1.97,2.53,2.83,3.32,3.7,4.09]

        
def getExperimentalKilns():
    return [
        BarrKilnT4(),
        BarrKilnT8(),
        BarrKilnT2(),
        BarrKilnT3(),
        BarrKilnT5(),
        BarrKilnT6(),
        BarrKilnT7(),
        BarrKilnT1(),
        BarrKilnT9(),
        TschengKilnA11(),
        TschengKilnA12(),
        TschengKilnA13(),
        TschengKilnA14(),
        TschengKilnA15(),
        TschengKilnA16(),
        TschengKilnA17(),
        TschengKilnA18(),
        TschengKilnA19(),
        TschengKilnA20(),
        TschengKilnA21(),
        TschengKilnA22(),
        TschengKilnA23(),
        TschengKilnA24(),
        TschengKilnA25(),
        TschengKilnA26(),
        TschengKilnA27(),
        TschengKilnA28(),
        TschengKilnA29(),
        TschengKilnA30(),
        TschengKilnA31(),
        TschengKilnA32(),
        TschengKilnA33(),
        TschengKilnA34(),
        TschengKilnA35(),
        TschengKilnA36(),
        TschengKilnA37(),
        TschengKilnA38(),
        TschengKilnA39(),
        TschengKilnA40(),
        TschengKilnA41(),
        TschengKilnA42(),
        TschengKilnA43(),
        TschengKilnA44(),
        TschengKilnA45(),
        TschengKilnA46(),
        TschengKilnA47(),
        TschengKilnA48(),
        TschengKilnA49(),
        TschengKilnA50(),
        TschengKilnA51(),
        TschengKilnA52(),
        TschengKilnA53(),
        TschengKilnA54(),
        KDOKilnTP21(), KDOKilnTP21r(),
        #KDOKilnTP26(),
    ]
