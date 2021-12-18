def load_defaults():
    print('Loading defaults...')
    # Framework
    days = 365 * 3 # Two years
    dt   = 0.01 # units: days    
    #url = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst2Agg.htmlTable?sst%5B(2020-04-26T12:00:00Z):1:(2020-04-26T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(47.135):1:(46.573)%5D%5B(-2.156):1:(-2.753)%5D,anom%5B(2020-04-26T12:00:00Z):1:(2020-04-26T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(47.135):1:(46.573)%5D%5B(-2.156):1:(-2.753)%5D,err%5B(2020-04-26T12:00:00Z):1:(2020-04-26T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(47.135):1:(46.573)%5D%5B(-2.156):1:(-2.753)%5D,ice%5B(2020-04-26T12:00:00Z):1:(2020-04-26T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(47.135):1:(46.573)%5D%5B(-2.156):1:(-2.753)%5D'
    
    # Parameters
    par = {}
    par['xk'] = 0.5
    par['Ingmax'] = 22 # units = mg POM h−1 g−1 
    par['ae'] = 0.25
    par['b'] = 0.66
    par['br'] = 0.8
    par['Ing1'] = 0.06
    par['Ing2'] = 19
    par['Shellcoef'] = 0.05 
    par['mgO2pom'] = 0.7
    par['Rgs1'] = 0.2
    par['sm'] = 1.5
    par['resp1'] = 0.2
    par['resp2'] = 0.1
    par['cms1'] = 0.0293 #maximal somatic growth parameter 1
    par['cms2'] = 0.044 #Maximal somatic growth parameter 2
    par['Shellgam'] = 0.5
    
    
    

    
    # Initial conditions
    InitCond = {}
    InitCond['SOMA'] = 0.3
    InitCond['RESGON'] = 0.06
    InitCond['SHELL'] = 0.08


    return  days, dt, par, InitCond
    
def run(days, dt, par, InitCond):
    print('Running model...')
    # Import libraries
    import numpy as np
    
    # Setup the framework 
    NoSTEPS = int(days / dt) # Calculates the number of steps 
    time = np.linspace(0,days,NoSTEPS) # Makes vector array of equally spaced numbers 
    
    # Create arrays of zeros
    B = np.zeros((NoSTEPS,),float) # Biomass 
    SOMA = np.zeros((NoSTEPS,),float) 
    RESGON = np.zeros((NoSTEPS,),float)
    SHELL = np.zeros((NoSTEPS,),float)
    
    
    # Initializing with initial conditions
    SHELL[0] = InitCond['SHELL']
    SOMA[0] = InitCond['SOMA'] 
    RESGON[0] = InitCond['RESGON'] 
    B[0] = InitCond['SOMA'] + InitCond['RESGON'] + InitCond['SHELL']
    
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
       
        #Forcing
        TEMP = 20 #forced average temperature, analytic
        X = 0.5   #Food proxy, forced analytic
       
        #SP = 0 #Fix later
        #DTM = 0 #Fix later
        f = X / (X + par['xk']) #Functional response, 
        INGtemp = par['Ing1'] * ((TEMP - par['Ing2']) ** 2) # Effect of temperature on ingestion - Units: mg h−1 g−1
        
        INGorg = (par['Ingmax'] - INGtemp) * f * (SOMA[t] ** par['b']) #Organic Ingestion - Units: mg h−1
        
        ABSorg = (par['ae'] * INGorg) / (1000.24) #Organic matter absorbtion - Units: Gd−1
        
       # if SOMA < 1.5 :
        #    maxSOMA = par['cms1'] #fix later - Units: Gd−1
        #else :
         #   maxSOMA = par['cms2'] 
            
        maxSOMA = par['cms2']
        
        RESP = (par['resp1'] ** (par['resp2'] * TEMP) * SOMA[t] ** par['br'] * par['mgO2pom']) / (1000 * 24) # Respiration - Units: G d−1
        
       
        
         # Shell gain - Units: Gd−1
        SHELLgain = (par['Shellcoef'] * ABSorg) * par['Shellgam']
        
        #SHELLgain constraint
        if ABSorg < RESP :
            SHELLgain = 0
            
        #Soma Gain - Units: Gd−1
        SOMAgain = min((maxSOMA), (ABSorg - SHELLgain)) - RESP
        
        # Reserve and gonad gain - Units: Gd−1
        RESGONgain = ABSorg - SHELLgain - SOMAgain 
        
        #
    
        
        # Spawning 
        #if RESGON[t]/B[t] < par['Rgs1']:
         #   Spawning = 0.
        #elif RESGON[t]/B[t] >= par['Rgs1']: #fix this
         #   Spawning = RESGON[t]
            #RESGON[t] = 0.
            
        #New spawning
        
      #  if RESGON/SOMA > par['Rgs1'] :
       #     Spawning = par['cs']
        
        # Update and step ------------------------------
        SOMA[t+1] = SOMA[t] + (SOMAgain * dt) 
        RESGON[t+1] = RESGON[t] + (RESGONgain * dt) #- Spawning
        SHELL[t+1] = SHELL[t] + (SHELLgain * dt)
        B[t+1] = SOMA[t+1] + RESGON[t+1] + SHELL[t+1]
        
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['INGorg'] = INGorg
    output['INGtemp'] = INGtemp
    output['ABSorg'] = ABSorg
    output['maxSOMA'] = maxSOMA
    output['time'] = time
    output['B'] = B 
    output['SOMA'] = SOMA
    output['RESGON'] = RESGON 
    output['SHELL'] = SHELL
    output['RESP'] = RESP
   
    
    
    print('Model run: DONE!!!')
    return  output

def plot(output):
    import matplotlib.pyplot as plt 
    # Plotting                      
    fig, (ax) = plt.subplots(1,1)
    ax.plot(output['time']/365,output['B'],'r-') 
    ax.plot(output['time']/365,output['SOMA'],'b-') 
    ax.plot(output['time']/365,output['RESGON'],'g.')
    ax.plot(output['time']/365,output['SHELL'], 'y-')
    #ax.plot(output['time']/365,output['INGorg'], 'm-')
    ax.legend(['B', 'SOMA', 'RESGON', 'SHELL']) 
    ax.set_ylabel('Mass (g)') 
    ax.set_xlabel('Time (years)')
    plt.show()          
    
    return

if __name__ == "__main__":
    print('Executing my_module.py')
    print('--------------------')
    
    days, dt, par, InitCond = load_defaults()
    output = run(days, dt, par, InitCond)
    plot(output)
    
    print('--------------------')