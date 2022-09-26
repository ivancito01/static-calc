                    ###################
                    #   STATIC CALC   #
                    ###################
#Version:0.0.1
#Autor: Morales,Ivan Ezequiel



import math
import numpy as np
import matplotlib.pyplot as plt


class MiembroArmadura():
    def __init__(self,elemento,area,mod_Elasticidad,nodo_Inicial,nodo_final,coordenada_Xi,coordenada_Yi,coordenada_Xf,coordenada_Yf,mom_Inercia,masa,aceleracion,tipo):
        
        ##################
        #DATOS DE ENTRADA#        
        ##################
        
        self.elem = elemento
        self.A = area #Area transversal de la sección
        self.E = mod_Elasticidad 
        self.xi = coordenada_Xi #coordenada en x del nodo inicial 
        self.yi = coordenada_Yi #coordenada en y del nodo inicial
        self.xf = coordenada_Xf #coordenada en x del nodo inicial 
        self.yf = coordenada_Yf #coordenada en y del nodo inicial
    