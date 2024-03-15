# NUM2.0

NUM 2.0 model by Grigoratou et al 2024. 
The model follows NUM 1.0's (Pompei et al., 2020) design, parameterazation and set up. The main differences compare to Pompei et al., 2020 can be summarized as follow:
NUM 2.0 does not include copepods' fecal pellet production
NUM 2.0 is more flexible compare to the NUM 1.0 on adding different plankton groups via a dataframe (e.g., the species dataframes in the folder).
NUM 2.0 represents the temperature effect on species uptake, clearence, ingestion rates via temperature optima following the Dutzkiewicz et al., (2015; 2020) Darwin- MIT approach.

The seasonal model can be run via: z_main_seasonal.m The ode solution can be found at: z_ode_copepod_model.m. Model's initial set up includes 200 species (112 protists, 64 active copepods and 24 passive copepods). 

For the Pompei et al., 2020 set up for number of plankton groups and sizes:
Choose “Pompei” for the “temp_option” in the z_ode_copepod_model.m
Choose “Pompei” for the “size_option” in the z_function_parameters.m
in the z_main_seasonal you need to change the following: 
nbrP=14; %number of protists size classes
nbrC_act=8; %number of copepod active feeders populations
nbrC_pass=3; %number of copepods passive feeders populations


Maria Grigoratou, March 2024
