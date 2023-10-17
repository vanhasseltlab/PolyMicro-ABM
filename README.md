# PolyMicro-ABM
## An Agent based modeling framework of an interactive, polymicrobial community during antimicrobial treatment

This folder contains all scripts and generated data needed to reproduce the analysis and figures shown in the manuscript ‘Pharmacodynamics of interspecies interactions’. All code was produced by Catharina Herzberg, LACDR, Leiden University, The Netherlands.

## Requirements/System:

The analysis has been performed with python 3.8. The python package mesa 1.1.1 was used to create the model (Model_v19.py). Further dependencies: matplotlib 3.5.3, numpy 1.23.3, pandas 1.4.4, seaborn 0.12.0


## Overview of files and order of execution:

The analysis was performed in 3 steps: 1. simulation, 2. data processing 3. plotting. The following lists all files in this software supplementary in the order of execution. 

1. Run 10 simulation scripts found in folder ‘Scripts/Simulation scripts’
- Input: Model_v19.py, Functions.py, PostSim_16.py
- Output in Data folder: 10 subfolders with the same name as scripts that generated them. Each contains db_agents.pkl, db_model.pkl, db_m_part.pkl, dict_model.pkl, settings.pkl and species.pkl

2A. Run Processing scripts found in folder ‘Scripts/Processing scripts’ in order: Create_Load_df_all_2.py, Create_Load_df_all_v17_2.py, Compare_rel_effect_log2.py
- Input: db_model of each Data subfolder
- Output in Data folder: df_all_v16&v17.pkl, df_all_s_mvnomv_v16&v17_VI.pkl, df_all_s_intrad1|0_v16&v17_VI.pkl

2B. Run Processing scripts found in folder ‘Scripts/Processing scripts’ in order: Create_df_dict_3_1B.py, Create_df_agents_part_1B_3.py
- Input: db_model of each of 8 Data subfolder with interactions present
- Output: df_dict_3.pkl in each of the 8 Data subfolder, df_agents_part_3.pkl, df_metrics_3.pkl

3A. Run plotting scripts found in folder ‘Scripts/Plotting scripts’ following 2A processing without single agent data
3A.1. PlotFiguresPart1_Improved.py
- Input:  df_all_v16&v17.pkl
- Output in Figures: Figure6_highdpi_v2.png, Figure7_highdpi_v2.png

3A.2. PlotFiguresPart2_Improved.py
- Input: df_all_v16&v17.pkl, df_all_s_mvnomv_v16&v17_VI.pkl
- Output: Figure8_highdpi.png, Figure9_highdpi_v5.png
	
3A.3. PlotFiguresPart3_Improved.py
- Input: df_all_s_intrad1|0_v16&v17_VI.pkl
- Output: Figure10-1_highdpi_v2.png, Figure10-2_highdpi_v2.png

3B. Run plotting scripts found in folder ‘Scripts/Plotting scripts’ following 2B processing with single agent data
3B. Plot_Figures_1B_Improved.py
- Input: df_metrics_3.pkl, df_agents_3.pkl
- Output: Figure4_highdpi.png, Figure3B_highdpi_allits.png, Figure3_highdpi_v4.png, Figure5_highdpi_v3.png
