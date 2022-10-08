                    ###################
                    #   STATIC CALC   #
                    ###################
#Version:0.0.1
#Autor: Morales,Ivan Ezequiel

import math
import numpy as np
import matplotlib.pyplot as plt


class MiembroArmadura():
    def __init__(self,elemento,area,mod_Elasticidad,nodo_Inicial,nodo_Final,coordenada_Xi,coordenada_Yi,coordenada_Xf,coordenada_Yf,mom_Inercia,masa,aceleracion,tipo):
        
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
        self.I = mom_Inercia
        self.M = masa
        self.G = aceleracion
        self.NI = nodo_Inicial
        self.NF = nodo_Final
        self.tipo = tipo
        #propiedades-resultados#

        self.L =self.Longitud()
        self.k_loc=self.Rig_Loc()
        self.Lx=self.Lambda_X()  #cos alpha#
        self.Ly=self.Lambda_Y()  #sen alpha#
        self.T=self.M_Transformacion()
        self.k_glob=self.Rig_Global()
        
        self.cosa = self.Lx
        self.sena = self.Ly
        self.q = ((self.M * self.G)/self.L)
        self.qx = self.q * self.sena
        self.qy = self.q * self.cosa
        
        
    def __str__(self):
        print("Area: ",self.A)
        print("Modulo Elastico:",self.E)
        print("Coordenada Inicial ({},{})".format(self.xi,self.yi))
        print("Coordenada Final ({},{})".format(self.xf,self.yf))
        print("Momenton de Inercia:",self.I)
        print("Masa:",self.M)
        print("Aceleracion:",self.G)
        print("Longitud:",self.L)
        print("Rigidez Local:\n",self.k_loc)
        print("Lx: ",self.Lx)
        print("Ly: ",self.Ly)
        print("Matriz De Transformacion\n",self.T)
        print("Matriz De rigidez global\n",self.k_glob)
        return ""
    
    #FUNCIONES DE LA CLASE#
    
    def Longitud(self):#calcula la longituda de la barra
        return math.sqrt((self.xf-self.xi)**2+(self.yf-self.yi)**2)

    def Rig_Loc(self):#calcula la matriz de rigidez de la barra       
        return np.array(
                        [[self.A*self.E/self.L,0,0,-self.A*self.E/self.L,0,0],
                        [0,12*(self.E*self.I/(self.L**3)),-6*(self.E*self.I/(self.L**2)),0,-(12*self.E*self.I/(self.L**3)),(-6*self.E*self.I/(self.L**2))],
                        [0,(-6*self.E*self.I/(self.L**2)),4*(self.E*self.I/(self.L)),0,6*(self.E*self.I/(self.L**2)),2*(self.E*self.I/(self.L))],
                        [-self.A*self.E/self.L,0,0,self.A*self.E/self.L,0,0],
                        [0,-12*(self.E*self.I/(self.L**3)),6*(self.E*self.I/(self.L**2)),0,12*(self.E*self.I/(self.L**3)),6*(self.E*self.I/(self.L**2))],
                        [0,-6*(self.E*self.I/(self.L**2)),2*(self.E*self.I/(self.L)),0,6*(self.E*self.I/(self.L**2)),4*(self.E*self.I/(self.L))]])
    
    #Lambda_X y Lambda_Y son los cosenos directores de las barras   
    def Lambda_X(self):
        return (self.xf-self.xi)/self.L

    def Lambda_Y(self):
        return (self.yf-self.yi)/self.L

        #############|Recordatorio|##########
        #self.Lx=self.Lambda_X()  #cos alpha#
        #self.Ly=self.Lambda_Y()  #sen alpha#
        #####################################

    def M_Transformacion(self):#Calcula la matriz para pasar de coordenadas localas a coordenadas globales
        return np.array(
            [[self.Lx,self.Ly,0,0,0,0],
             [-self.Ly,self.Lx,0,0,0,0],
             [0,0,1,0,0,0],
             [0,0,0,self.Lx,self.Ly,0],
             [0,0,0,-self.Ly,self.Lx,0],
             [0,0,0,0,0,1],])
        
    def Rig_Global(self):#pasaje de matriz de rigidez del elemento a coordenadas globales
        '''producto vectorial entre matriz t transpuesta y matriz k de rigidez local'''
        Ttrans_X_k_loc=np.matmul(np.transpose(self.T),self.k_loc)
        '''producto vectorial entre matriz t transpuesta con matriz k de rigidez local matriz t'''
        return np.matmul(Ttrans_X_k_loc,self.T)
    
    
print(MiembroArmadura("b1",1,100,1,2,0,0,7,5,1000,50,8.91,1))



class calculo_esfuerzos_Ex():
    def __init__(self,Ug,T,K_loc,Ni,Nf,fqe):
        self.Ug = Ug #desplazamientos en ejes globales
        self.T = T #matriz de transformacion del elemento
        self.K_loc = K_loc #matriz de transformacion en coordenadas locales
        self.ni = Ni #nodo inicial
        self.nf = Nf #nodo final
        self.fqe=fqe
        #funciones internas de la clase
        self.Ug_red = self.reduction_Vec_Ug()
        self.Ue = self.pasaje_Ejes_Loc()
        self.ext = self.esfuerzos_extremos()
        self.extg = self.esfuerzos_extremos_g()
    def reduction_Vec_Ug(self):
        vq=np.zeros((6,1))
        vq[0]=self.Ug[(self.ni*3-2)-1]
        vq[1]=self.Ug[(self.ni*3-1)-1]
        vq[2]=self.Ug[(self.ni*3)-1]
        vq[3]=self.Ug[(self.nf*3-2)-1]
        vq[4]=self.Ug[(self.nf*3-1)-1]
        vq[5]=self.Ug[(self.nf*3)-1]
        return vq
    def pasaje_Ejes_Loc(self):
        return np.matmul(self.T,self.Ug_red)
    def esfuerzos_extremos(self):#en ejes locales
        return np.add(np.matmul(self.K_loc,self.Ue),self.fqe)
    def esfuerzos_extremos_g(self):
        Trans=np.transpose(self.T)
        return np.matmul(Trans,self.ext)
        
class Vector_Fuerzas():

    def __init__(self,elemento,Q1,Q2,F,L,matDeTrans,masa,aceleracion,lx,ly):

      #Datos de entrada#
        self.elem = elemento        
        self.Q1 = Q1
        self.Q2 = Q2
        self.F  = F
        self.L  = L
        self.T  = matDeTrans
        self.M = masa
        self.G = aceleracion
        self.cosa = lx
        self.sena = ly
      #Propiedades - Resultados# 

        self.vFe = self.vector_cargas_distribuidas_loc()
        self.vFg = self.vector_cargas_distribuidas_glob()
        self.vFme = self.vector_carga_masica_distribida_loc()
        self.vFmg = self.vector_carga_masica_distribida_glob()
        self.vFte = self.vector_carga_total_loc()
        self.vFtg = self.vector_carga_total_glob()

    def __str__(self):
        print("vector de cargas distribuidas en ejes locales\n:",self.vFe)
        print("vector de cargas distribuidas en ejes globales\n:",self.vFg)
        print("vector de cargas masicas en ejes locales\n:",self.vFme)
        print("vector de cargas masicas en ejes globales\n:",self.vFmg)


        return""
    def vector_cargas_distribuidas_loc(self):        
        vq_loc =np.zeros((6, 1))
        q1=-self.Q1
        q2=-self.Q2
        F=self.F
        L=self.L
        if(abs(q1)<abs(q2)):
          vq_loc[0][0]  = F*L/2
          vq_loc[1][0]  = -((3/20)*(q2-q1)*L)-(q1*L/2)
          vq_loc[2][0]  = ((q2-q1)*(L**2)/30)+(q1*(L**2)/12)
          vq_loc[3][0]  = F*L/2
          vq_loc[4][0]  = -((7/20)*(q2-q1)*L)-(q1*L/2)
          vq_loc[5][0]  = -((q2-q1)*(L**2)/20)-(q1*(L**2)/12)
        elif(abs(q1)==abs(q2)):
          vq_loc[0][0]  = F*L/2
          vq_loc[1][0]  = -q1*L/2
          vq_loc[2][0]  = (q1)*(L**2)/12
          vq_loc[3][0]  = F*L/2
          vq_loc[4][0]  = -(q1*L/2)
          vq_loc[5][0]  = -(q1)*(L**2)/12
        elif(abs(q1)>abs(q2)):
          vq_loc[0][0]  = F*L/2
          vq_loc[1][0]  = -((3/20)*(q1-q2)*L)-(q2*L/2)
          vq_loc[2][0]  = ((q1-q2)*(L**2)/30)+(q2*(L**2)/12)
          vq_loc[3][0]  = F*L/2
          vq_loc[4][0]  = -((7/20)*(q1-q2)*L)-(q2*L/2)
          vq_loc[5][0]  = -((q1-q2)*(L**2)/20)-(q2*(L**2)/12)        
        return vq_loc
    
    def vector_cargas_distribuidas_glob(self):
        Trans=np.transpose(self.T)
        vector2=np.matmul(Trans,self.vFe)
        return vector2
    def vector_carga_masica_distribida_loc(self):
        q = ((self.M * self.G)/self.L)
        qx = -q * self.sena
        qy = -q * self.cosa
        vq_loc =np.zeros((6, 1))
        vq_loc[0][0]  = -(qx*self.L)/2
        vq_loc[1][0]  = -(qy*self.L)/2
        vq_loc[2][0]  = ((qy)*(self.L**2))/12
        vq_loc[3][0]  = -(qx*self.L)/2
        vq_loc[4][0]  = -(qy*self.L/2)
        vq_loc[5][0]  = -(qy)*(self.L**2)/12
        return vq_loc
    def vector_carga_masica_distribida_glob(self):
        Trans=np.transpose(self.T)
        vector2=np.matmul(Trans,self.vFme)
        return vector2
    def vector_carga_total_loc(self):
      return np.add(self.vFe,self.vFme)

    def vector_carga_total_glob(self):
       return np.add(self.vFg,self.vFmg)
        
class AnalisisMatricial():
   
    def __init__(self, tbl_Elem, tbl_Nods, tbl_Frza,tbl_Rest,tbl_Q_distr):
        
        ####################
        # Datos de entrada #
        ####################
        
        self.tE = tbl_Elem # Tabla de elementos.
        self.tN = tbl_Nods # Tabla de nodos.
        self.tF = tbl_Frza # Tabla de fuerzas.
        self.tR = tbl_Rest # Tabla de desplazamientos.
        self.tQ = tbl_Q_distr# Tabla de cargas distribuidas.
        
        
        ###############
        # Propiedades #
        ###############
        
        self.nE = len(self.tE) # Número de elementos en la estructura.
        self.nN = len(self.tN) # Número de nodos en la estructura.
        self.nF = len(self.tF)# Número de fuerzas puntuales.   
        self.nQ = len(self.tQ)# Número de cargas deistribuidas.     
        self.nGl = self.nN*3 # Número de grados de libertad en la estructura.      
        self.N = self.Diccionario_de_nodos() # Contiene posision de cada nodo y sus grados de libertas
        self.Dtipo = self. tipo_elemento() # Actualmente esta funcion no es utilizada(actualizacion pendiente)
        self.Armad = self.Armadura() # Contiene la informacion mas relevante de cada elemento
        self.mgdl = self.matriz_GL() # Es una matriz de filas=elemento columnas=grados de livertad del elemento
        self.mkg = self.matriz_K_global() # Matriz de rigidez global de toda la estructura
        self.vq = self.vector_cargas_puntuales()  # Vector de esfuerzos puntuales en nodos
        self.vqd = self.Vectores_fuerzas() # Vector de esfuerzos nodales producidos por cargas distribuidas 
        self.vqd_loc = self.Vectores_fuerzas_loc()
        self.VFG = self.matriz_F_global() # Vector esfuerzos producidos por las cargas de toda la estructura
        self.vGDE = self.grados_libertad_eliminados() # Grados de libertad asociados a las reacciones
        self.mkgr = self.matriz_reducida() # Reduccion de la matriz k glob (eliminacion de filas y columnas asociadas a gdl de reacciones)
        self.mkgri = self.matriz_reducida_inversa() # Inversion de la matriz k glob reducida
        self.Fqg =  self.vec_cargas_reducida() # Reduccion del vector de fuerzas
        self.Fg = self.vec_fuerzas_reducida() # Reduccion del vector de fuerzas
        self.DesD= self.desplazamientos_desconocidos() # Resolvemos para encontrar los desplazamientos
        self.vecDes =self.des_ampliados() # Rearmamos el vector desplazamiento con los dezplazamientos que conocemos
        self.KgUg =self.Kg_X_Des_ampl() #producto de la matriz k glob con los desplazamientos en coordenadas globales
        self.veds = self.vector_de_esfuerzos_desc()
        self.extrem=self.esfuerzos_extremos()

    
              
    def __str__(self):
#         print("Numero de Elementos: ", self.nE)
#         print("Numero de Nodos: ", self.nN)
#         print("Numero de Grados de Libertad: ", self.nGl)
#         #print("matriz de grados de libertad:\n",self.mgdl)        
#         print("diccionario de nodos:\n",self.N)    
#         #print("elemento [3]:\n\n",self.Armad[2].k_glob,"\n",self.Armad[2].k_loc)        
#         for i in range(self.nE):
#           print("vector de cargas del elemto [",i+1,"]:\n",self.vqd[1])
#         print("vector de fuerzas global:\n",self.VFG)
#         print("matriz de rigidez global:\n",self.mkg)
#         print("vector de grados de libertad eliminados:\n",self.vGDE)
#         print("vector de cargas:\n",self.Fqg) 
#         print("vector de fuerzas:\n",self.Fg)
#         print("vector de desplazamientos:\n",self.vecDes)    
#         print("vector de fuerzas netas:\n",self.KgUg)      
        print("\n",self.extrem[0].ext)
        print("reducido\n",self.extrem[0].Ug_red)
        print("local\n",self.extrem[2].Ue)
        return ""
    
    def Diccionario_de_nodos(self):
      DiccionarioNodos = {}
      #NODO_key, [cx, cy, glx, gly, glz]
      for i in range(self.nN):            
            NODO_key = self.tN[i][0]
            cx = self.tN[i][1]
            cy = self.tN[i][2]
            tipo = self.tN[i][3]
            DiccionarioNodos.setdefault(NODO_key, [cx, cy,(i+1)*3-2,(i+1)*3-1,(i+1)*3,tipo])

      return DiccionarioNodos
    def tipo_elemento(self):
        tipo = []
        for i in range(self.nE):
            ni=self.tE[i][3]
            nf=self.tE[i][4]
            tni=self.N[ni][5]
            tnf=self.N[nf][5]
            if(tni=="Rigido" and tnf=="Rigido"):
              tipo.append(1)
            elif(tni=="Articulado" and tnf=="Articulado"):
              tipo.append(2)
            elif(tni=="Articulado" and tnf=="Rigido"):
              tipo.append(3)
            elif(tni=="Rigido" and tnf=="Articulado"):
              tipo.append(4)
        return tipo

    def Armadura(self):
        #Calcula la lista de propiedades de cada elemento de la armadura.#
        Elem = []
        for i in range(self.nE):                
            el = self.tE[i][0]
            a = self.tE[i][1]
            me = self.tE[i][2]
            NI = self.tE[i][3]
            NF = self.tE[i][4]
            xi = self.N[NI][0]
            yi = self.N[NI][1]
            xf = self.N[NF][0]
            yf = self.N[NF][1]             
            mom = self.tE[i][5]
            M = self.tE[i][6]
            G = self.tE[i][7]
            tip = self.Dtipo[i]
            Elem.append(MiembroArmadura(el, a, me,NI,NF, xi, yi, xf, yf,mom,M,G,tip))
        return Elem

    def Vectores_fuerzas(self):
        #Calcula la lista de vectores de fuerzas de cada elemento de la armadura.#
        vec_fuerzas = []
        for i in range(self.nE):
          masa = self.Armad[i].M
          lx = self.Armad[i].Lx
          ly = self.Armad[i].Ly
          L = self.Armad[i].L
          T = self.Armad[i].T
          aceleracion = self.Armad[i].G
          vec_fuerzas.append(Vector_Fuerzas(i,0,0,0,L,T,masa,aceleracion,lx,ly).vFmg)

        for i in range(self.nQ):
            elemento = int(self.tQ[i][0])-1                
            Q1  = self.tQ[i][1]
            Q2  = self.tQ[i][2]
            F = self.tQ[i][3]
            L = self.Armad[elemento].L
            T = self.Armad[elemento].T
            masa = self.Armad[elemento].M
            aceleracion = self.Armad[elemento].G  
            lx = self.Armad[elemento].Lx
            ly = self.Armad[elemento].Ly                             
            vec_fuerzas[elemento]=Vector_Fuerzas(elemento,Q1,Q2,F,L,T,masa,aceleracion,lx,ly).vFtg
        return vec_fuerzas
    def Vectores_fuerzas_loc(self):
        #Calcula la lista de vectores de fuerzas de cada elemento de la armadura.#
        vec_fuerzas = []
        for i in range(self.nE):
          masa = self.Armad[i].M
          lx = self.Armad[i].Lx
          ly = self.Armad[i].Ly
          L = self.Armad[i].L
          T = self.Armad[i].T
          aceleracion = self.Armad[i].G
          vec_fuerzas.append(Vector_Fuerzas(i,0,0,0,L,T,masa,aceleracion,lx,ly).vFme)

        for i in range(self.nQ):
            elemento = int(self.tQ[i][0])-1                
            Q1  = self.tQ[i][1]
            Q2  = self.tQ[i][2]
            F = self.tQ[i][3]
            L = self.Armad[elemento].L
            T = self.Armad[elemento].T
            masa = self.Armad[elemento].M
            aceleracion = self.Armad[elemento].G  
            lx = self.Armad[elemento].Lx
            ly = self.Armad[elemento].Ly                             
            vec_fuerzas[elemento]=Vector_Fuerzas(elemento,Q1,Q2,F,L,T,masa,aceleracion,lx,ly).vFte
        return vec_fuerzas
        
    def matriz_GL(self):
        M=[]
        for i in range(self.nE):
          M.append([0]*6)
        for i in range(self.nE):
          Dn=self.tE[i][3] #nodo inicial
          An=self.tE[i][4] #nodo final
          M[i][0]=(Dn*3)-2
          M[i][1]=(Dn*3)-1
          M[i][2]=(Dn*3)
          M[i][3]=(An*3)-2
          M[i][4]=(An*3)-1
          M[i][5]=(An*3)
       
        return np.array(M)
    def matriz_K_global(self):      
        K = np.zeros((self.nGl, self.nGl))
        for e in range(self.nE):
            ke_global = self.Armad[e].k_glob
            for i in range(6):
             for j in range(6):
              a = self.mgdl[e,i]-1
              b = self.mgdl[e,j]-1
              K[a,b] = ke_global[i,j] + K[a,b]
        return K
        
    def vector_cargas_puntuales(self):
        vq = np.zeros((self.nGl, 1))
        for i in range(self.nF):
          grado=0
          Force=self.tF[i][0]  #magnitud de la fuerza
          Nodo=self.tF[i][1]  #nodo donde se aplica la fuerza          
          Direction=self.tF[i][2]  #direccion de aplicacion de la fuerza
          #calculo del grado de libertad de aplicacion de la fuerza#
          if(Direction=="DX"):#Direccion de X#
              grado=(Nodo*3)-2
          elif(Direction=="DY"):#Direccion de Y#
              grado=(Nodo*3)-1
          elif(Direction=="MA"):#Momento Aplicado#
              grado=(Nodo*3)
          vq[grado-1][0]=Force
        return vq

    def matriz_F_global(self):
        vq = np.zeros((self.nGl, 1))
        for i in range(self.nE): 

          ni=self.Armad[i].NI
          nf=self.Armad[i].NF        
          
          vq[int((ni*3)-2)-1] = vq[int((ni*3)-2)-1] + self.vqd[i][0]          
          vq[int((ni*3)-1)-1] = vq[int((ni*3)-1)-1] + self.vqd[i][1]
          vq[int((ni*3))-1]   = vq[int((ni*3))-1] + self.vqd[i][2]
          vq[int((nf*3)-2)-1] = vq[int((nf*3)-2)-1] + self.vqd[i][3]         
          vq[int((nf*3)-1)-1] = vq[int((nf*3)-1)-1] + self.vqd[i][4]
          vq[int((nf*3))-1]   = vq[int((nf*3))-1] + self.vqd[i][5]      
                       
        return vq      
    
    def grados_libertad_eliminados(self):
      vec=[]
      for i in range(len(self.tR)):
        nodo = int(self.tR[i][0])
        Rx = self.tR[i][1]
        Ry = self.tR[i][2]
        Rm = self.tR[i][3]
        if(Rx=="Si"):
          vec.append((nodo*3-2)-1)
        if(Ry=="Si"):
          vec.append((nodo*3-1)-1)
        if(Rm=="Si"):
          vec.append((nodo*3)-1)
      return vec
    def matriz_reducida(self):
      mat=self.mkg
      mat1=np.delete(mat,(self.vGDE), axis = 0)
      mat2=np.delete(mat1,(self.vGDE), axis = 1)
      return mat2
    def matriz_reducida_inversa(self):
      return np.linalg.inv(self.mkgr)
    def vec_cargas_reducida(self):
      return np.delete(self.VFG,(self.vGDE), axis = 0)
    def vec_fuerzas_reducida(self):
      return np.delete(self.vq,(self.vGDE), axis = 0)
    def desplazamientos_desconocidos(self):
       vec_red = self.Fg-self.Fqg
       return np.matmul(self.mkgri,vec_red)
    def des_ampliados(self):
      vec=self.DesD
      for i in self.vGDE:
        vec=np.insert(vec,i,0,axis = 0)
      return vec
    def Kg_X_Des_ampl(self):
       return np.matmul(self.mkg,self.vecDes)
    def vector_de_esfuerzos_desc(self):       
        return np.add(self.KgUg,self.VFG)

    def esfuerzos_extremos(self):
        vq=[]
        for i in range(self.nE):
            Ug=self.vecDes
            T=self.Armad[i].T
            K_loc=self.Armad[i].k_loc
            Ni=self.Armad[i].NI
            Nf=self.Armad[i].NF
            fqe=self.vqd_loc[i]
            vq.append(calculo_esfuerzos_Ex(Ug,T,K_loc,Ni,Nf,fqe))
        return vq
    
    
#################################
#P R U E B A DEL P R O G R A M A#
#################################




A1 = 0.0003738
A2 = 0.0004378 
A3 = 0.001536
A4 = 0.000675
A5 = 0.000444
A6 = 0.000624
A7 = 0.001954
A8 = 0.02379
A9 = 0.0004378
A10= 0.000784
A11= 0.000775
A12= 0.000875
A13= 0.0011299
A14= 0.000975
Et = 210000000000 # Modulo de Elasticidad en, kgf/m2
Ec = 21000000000
I1 = 0.00000031751#80x60x1.6
I2 = 0.00000041587#80x40x1.6
I3 = 0.000004040990#150x50x4
I4 = 0.0000006264100 #80x60x2.5
I5 = 0.000000511450 #80x60x2
I6 = 0.0000010773100 #120x40x2
I7 = 0.0000038521900 #140x60x5.51
I8 = 0.0000060172600 #150x50x6.35
I9 = 0.0000004870000 #90x50x1.6
I10= 0.0000021326500 #150x50x2
I11= 0.0000010703600 #100x60x2.5
I12= 0.0000016693200 #100x50x2.5
I13= 0.0000008681000 #40x80x5.15
I14= 0.0000024432000 #140 x60x2.5
G = (2.5 * 9.8)

# Tabla de Elementos:
tbl_Elem = [
    ['E1', A10, Ec,1,2,I10,72.7,G],
    ['E2', A3, Et,2,3,I3,78.59,G],
    ['E3', A10, Et,3,4,I10,72.88,G],
    ['E4', A10, Et,4,5,I10,66.29,G],
    ['E5', A4, Ec,5,6,I4,66,G],
    ['E6', A8, Ec,6,7,0.0000009687300,25.1,G],
    ['E7', A8, Ec,1,8,I8,28.5,G],
    ['E8', A8, Ec,2,9,I8,33,G],
    ['E9', A2, Ec,9,3,I2,21,G],
    ['E10', A7, Ec,3,10,I7,28,G],
    ['E11', 0.00855, Et,4,14,0.0000007791900,18.33,G],
    ['E12', A8, Et,5,11,I8,33.33,G],
    ['E13', A8, Ec,6,11,0.0000009687300,25.34,G],
    ['E14', A8, Ec,7,12,I8,18.33,G],
    ['E15', A2,Et,8,9,I2,15.4,G],
    ['E16', A3, Ec,9,10,I3,18.4,G],
    ['E17', A8, Et,11,12,0.0000009687300,19.50786,G],
    ['E18', A11, Ec,10,13,I11,20,G],
    ['E19', A9, Ec,15,11,I9,12,G],
    ['E20', A9, Ec,13,14,I9,18,G],
    ['E21', A9, Ec,14,15,I9,12,G]]
# Tabla de Nodos:
#["id",x,y,Restriccion en X,Restriccion en Y, Restriccion de momemto]
#las restricciones se introducen como "Si" y "No"
# Tabla de Nodos:
tbl_Nods = [
    [1, 0, 0, "Rigido"],
    [2, 0.75, 0, "Rigido"],
    [3, 1.714, 0, "Rigido"],
    [4, 2.785,0, "Rigido"],
    [5, 3.261, 0, "Rigido"],
    [6, 3.820, 0, "Rigido"],
    [7, 4.570, 0.457, "Rigido"],
    [8, 0, 0.476, "Rigido"],
    [9, 1.335, 0.526, "Rigido"],
    [10, 1.714, 0.626, "Rigido"],
    [11, 4.237, 0.814, "Rigido"],
    [12, 4.570, 0.764, "Rigido"],
    [13, 2.216, 1.290, "Rigido"],
    [14, 2.931, 1.310, "Rigido"],
    [15, 3.644, 1.290, "Rigido"]]
# Tabla de Fuerzas:
tbl_Frza = [
    #[ -20, 7, 'DY']
    ]
# Tabla de Cargas Distribuidas
#[Elemento,Q1,Q2]
tbl_Q_distr =[
    #  [1,0,5.7143,0],
    #  [2,5.7143,10,0],
    #  [3,8,8,0],
    #  [4,8,8,0],
    #  [7,8,8,0], 
    #  [8,8,8,0],           
              ]


# Tabla de Desplazamientos:
tbl_Desp = [
            [2, 'Si', 'Si', 'No'],
            [6, 'Si', 'Si', 'No']]
AE = AnalisisMatricial(tbl_Elem,tbl_Nods, tbl_Frza, tbl_Desp,tbl_Q_distr)
All =AnalisisMatricial(tbl_Elem,tbl_Nods, tbl_Frza, tbl_Desp,tbl_Q_distr).Armad
print("Carga Existosa")

