import numpy as np

class Gas:

    def __init__(self,
                CH4 = 0,
                C2H6 = 0,
                C3H8 = 0,
                C4H10 = 0,
                C5H12 = 0,
                N2 = 0,
                CO2 = 0,
                O2 = 0,
                CO = 0,
                H2 = 0,
                alpha = 1):
        self.CH4 = CH4
        self.C2H6 = C2H6
        self.C3H8 = C3H8
        self.C4H10 = C4H10
        self.C5H12 = C5H12
        self.N2 = N2
        self.CO2 = CO2
        self.O2 = O2
        self.CO = CO
        self.H2 = H2
        self.alpha = alpha

    @staticmethod
    def m_to_n(x):
        return (2 * x + 2)

    @property
    def heat_of_combustion(self):
        combustion_products = np.array([self.CH4,self.C2H6,self.C3H8,self.C4H10])  
        coeff = np.array([358.2,637.46,860.05,1185.8])  
        return sum(coeff*combustion_products)

    @property
    def air_volume(self):
      combustion_products = np.array([self.CH4,self.C2H6,self.C3H8,self.C4H10,self.C5H12]) 
      V0 = []
      for indx, value in enumerate(combustion_products):
          m = np.array(indx + 1)
          n = np.array(self.m_to_n(m))
          V0.append(m+n/4 * value)
      return 0.0476 * sum(V0)

    @property
    def triatomic_gases_volume (self):
       combustion_products = [self.CH4,self.C2H6,self.C3H8,self.C4H10,self.C5H12]
       V0RO2 = []
       for indx, value in enumerate(combustion_products):
           m = np.array(indx + 1)
           n = np.array(self.m_to_n(m))
           V0RO2.append(m*value)
       V0RO2 = 0.01*(self.CO2 + sum(V0RO2))
       return V0RO2

    @property
    def water_volume (self):
       combustion_products = [self.CH4,self.C2H6,self.C3H8,self.C4H10,self.C5H12]
       V0HO2 = []
       for indx, value in enumerate(combustion_products):
           m = np.array(indx + 1)
           n = np.array(self.m_to_n(m))
           V0HO2.append((n/2)*value)
       V0HO2 = 0.01*(1.61*self.air_volume + sum(V0HO2))
       return V0HO2

    @property
    def nitrogen_volume(self):
        return 0.79 * self.air_volume + 0.01 * self.N2
    
    @property
    def actual_volume_of_water_vapor(self):
        return self.water_volume + 0.0161 * (self.alpha - 1)

    @property
    def full_volume_of_combustion_products(self):
        return self.triatomic_gases_volume + self.nitrogen_volume + self.actual_volume_of_water_vapor +(self.alpha - 1) * self.air_volume

    
    def heat_capacity_CO2(self,t):
        temp = np.array([t ** 3, t ** 2, t, 1])
        coeff = np.array([4.5784 * 10 ** (-11), -1.51719 * 10 ** (-7), 2.50113 * 10 ** (-4), 0.382325])
        return 4.1868*(sum(temp * coeff))

    
    def heat_capacity_N2(self,t):
        temp = np.array([t ** 3, t ** 2, t, 1])
        coeff = np.array([-2.24553 * 10 ** (-11),4.85082 * 10 ** (-8),-2.90598 * 10 ** (-6),0.309241])
        return 4.1868*(sum(temp*coeff))
 
    
    def heat_capacity_H20(self,t):
        temp = np.array([t ** 3, t ** 2, t, 1])
        coeff = np.array([-2.10956 * 10 ** (-11), 4.9732 * 10 ** (-8), 2.60629 * 10 ** (-5), 0.356691])
        return 4.1868*(sum(temp*coeff))

    
    def heat_capacity_air(self,t):
        temp = np.array([t ** 3, t ** 2, t, 1])
        coeff = np.array([-2.1717 * 10 ** (-11), 4.19344 * 10 ** (-8), 8.00891 * 10 ** (-6), 0.315027])
        return 4.1868*(sum(temp*coeff))

 
    def enthalpy_of_pure_combustion_products(self,t):
        return (self.triatomic_gases_volume * self.heat_capacity_CO2(t) + self.nitrogen_volume * self.heat_capacity_N2(t) + self.water_volume * self.heat_capacity_H20(t)) * t

   
    def enthalpy_air(self,t):
        return (self.stoichiometric_air_expenditure * self.heat_capacity_air(t=t))

    def relative_enthalpy(self,t):
        return self.enthalpy_of_pure_combustion_products(t) + (self.alpha - 1) / self.enthalpy_air(t)

    def find_PV_RT (p = None,
                    v = None,
                    t = None):
        R = 8.314 
        if p == None:
            return (R*t)/v
        elif v == None:
            return (R*t)/p
        elif t == None:
            return (p*v)/R










