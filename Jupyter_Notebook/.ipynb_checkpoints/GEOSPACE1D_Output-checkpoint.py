# -*- coding: utf-8 -*-
"""
Created on 10/29/2019
Last modified on 19/12/2022

@author: Concetta D'Amato
@license: creative commons 4.0
"""


import os
import warnings
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.style as style
import matplotlib.dates as mdates
import plotly.express as px
import plotly.graph_objects as go
from IPython.display import Image

from timeseries_writer import *

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')

def show_var(variabile):
    df = pd.read_csv(variabile+'_1.csv',skiprows=6, sep=',', parse_dates=[0], na_values=-9999,usecols=[1,2])
    df.columns = ['Datetime',variabile]

    fig = px.line(df, x='Datetime', y=variabile)
    if variabile=='airT':
        fig.update_yaxes(title_text='Temperature [°C]')
    if variabile=='Wind':
        fig.update_yaxes(title_text='Velocità vento [m/s] ')
    if variabile=='GHF':
        fig.update_yaxes(title_text='Flussio di calore dal suolo')
    if variabile=='Pres':
        fig.update_yaxes(title_text='Pressione')
    if variabile=='SoilMoisture_sin':
        fig.update_yaxes(title_text='Contenuto di acqua')
    if variabile=='LAI':
        fig.update_yaxes(title_text='LAI')
    if variabile=='RH':
        fig.update_yaxes(title_text='Umidità relativa')
    if variabile=='ShortwaveDirect':
        fig.update_layout(title='Short wave direct')
        fig.update_yaxes(title_text='Radiazione corta diretta [$W \cdot m^{−2}$]')
    if variabile=='ShortwaveDiffuse':
        fig.update_layout(title='Short wave diffuse')
        fig.update_yaxes(title_text='Radiazione corta diffusa [$W \cdot m^{−2}$]')
    if variabile=='LongDownwelling':
        fig.update_layout(title='Long wave')
        fig.update_yaxes(title_text='Radiazione lunga incidente [$W \cdot m^{−2}$]')
    if variabile=='Net':
        fig.update_layout(title='Net radiation')
        fig.update_yaxes(title_text='Radiazione netta [$W \cdot m^{−2}$]')
    
    #fig.update_xaxes(rangeslider_visible=True)
    fig.show()

    
def show_stress(a,b,c,d):
    warnings.filterwarnings('ignore')
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1) 
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'Potential'  
    kl.Potential[kl.Potential<0]=0
      
    kl3 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl3 = kl3.drop(['Format'],axis=1)
    kl3.columns.values[0] = 'Date'
    kl3.columns.values[1] = 'Environmental_Stress'
    kl3.Environmental_Stress[kl3.Environmental_Stress<0]=0
        
    kl2 = pd.read_csv(c,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'Water_Stress'
    kl2.Water_Stress[kl2.Water_Stress<0]=0
    
    kl4 = pd.read_csv(d,skiprows=6,parse_dates=[1])
    kl4 = kl4.drop(['Format'],axis=1)
    kl4.columns.values[0] = 'Date'
    kl4.columns.values[1] = 'Total_Stress'
    kl4.Total_Stress[kl4.Total_Stress<0]=0
    
    
   
    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['Potential'], mode='lines', name='Potential'))
    fig.add_trace(go.Scatter(x=kl3['Date'], y=kl3['Environmental_Stress'], mode='lines', name='Environmental_Stress'))
    fig.add_trace(go.Scatter(x=kl2['Date'], y=kl2['Water_Stress'], mode='lines', name='Water_Stress'))
    fig.add_trace(go.Scatter(x=kl4['Date'], y=kl4['Total_Stress'], mode='lines', name='Total_Stress'))
    

    fig.update_xaxes(rangeslider_visible=True)
    
    fig.update_layout(
        #title= a,
        #xaxis_title="Date",
        #yaxis_title= a,
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        legend_title="Active Stress",
        font=dict(size=14))
    
    fig.show()
    
    
def show_stress2(a,b):
    warnings.filterwarnings('ignore')
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1) 
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'Potential'  
    kl.Potential[kl.Potential<0]=0
        
    kl2 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'Water_Stress'
    kl2.Water_Stress[kl2.Water_Stress<0]=0

    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['Potential'], mode='lines', name='Potential'))
    fig.add_trace(go.Scatter(x=kl2['Date'], y=kl2['Water_Stress'], mode='lines', name='Water_Stress'))
    

    fig.update_xaxes(rangeslider_visible=True)
    
    fig.update_layout(
        #title= a,
        #xaxis_title="Date",
        #yaxis_title= a,
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        legend_title="Active Stress",
        font=dict(size=14))
    
    fig.show()



def show_E_T(a,b,c):
    warnings.filterwarnings('ignore')
    ############################ GRAFICO 1 ########################################
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1)
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'EvapoTranspiration'
    kl.EvapoTranspiration[kl.EvapoTranspiration<0]=0
   
    
    kl2 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'Evaporation'
    kl2.Evaporation[kl2.Evaporation<0]=0
  
    
    kl3 = pd.read_csv(c,skiprows=6,parse_dates=[1])
    kl3 = kl3.drop(['Format'],axis=1)
    kl3.columns.values[0] = 'Date'
    kl3.columns.values[1] = 'Transpiration'
    kl3.Transpiration[kl3.Transpiration<0]=0
    
    
    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['EvapoTranspiration'], mode='lines', name='EvapoTranspiration'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl2['Evaporation'], mode='lines', name='Evaporation'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl3['Transpiration'], mode='lines', name='Transpiration'))

   #fig.update_xaxes(rangeslider_visible=True)
    
    
    fig.update_layout(
        title='Compare Evaporation and Traspiration fluxes',
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        #xaxis_title="Date",
        yaxis_title="[$W\cdot m^{−2}$]",
        #legend_title="Date",
        font=dict(size=14))
    fig.show()
    
    
def show_E_T_(a,b,c,d):
    warnings.filterwarnings('ignore')
    ############################ GRAFICO 1 ########################################
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1)
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'EvapoTranspiration'
    kl.EvapoTranspiration[kl.EvapoTranspiration<0]=0
   
    
    kl2 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'Evaporation'
    kl2.Evaporation[kl2.Evaporation<0]=0
  
    
    kl3 = pd.read_csv(c,skiprows=6,parse_dates=[1])
    kl3 = kl3.drop(['Format'],axis=1)
    kl3.columns.values[0] = 'Date'
    kl3.columns.values[1] = 'Transpiration'
    kl3.Transpiration[kl3.Transpiration<0]=0
    
    
    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['EvapoTranspiration'], mode='lines', name='EvapoTranspiration'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl2['Evaporation'], mode='lines', name='Evaporation'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl3['Transpiration'], mode='lines', name='Transpiration'))

   #fig.update_xaxes(rangeslider_visible=True)
    
    
    fig.update_layout(
        title='Compare Evaporation and Traspiration fluxes',
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        #xaxis_title="Date",
        yaxis_title=d,
        #legend_title="Date",
        font=dict(size=14))
    fig.show()
    
def show_compare_EvapoTranspiration(a,b,c):
    warnings.filterwarnings('ignore')
    ############################ GRAFICO 1 ########################################
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1)
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'Prospero_ET'
    kl.Prospero_ET[kl.Prospero_ET<0]=0
  
    kl2 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'Pristley_Taylor_ET'
    kl2.Pristley_Taylor_ET[kl2.Pristley_Taylor_ET<0]=0
  
    kl3 = pd.read_csv(c,skiprows=6,parse_dates=[1])
    kl3 = kl3.drop(['Format'],axis=1)
    kl3.columns.values[0] = 'Date'
    kl3.columns.values[1] = 'PenmanMontheithFAO_ET'
    kl3.PenmanMontheithFAO_ET[kl3.PenmanMontheithFAO_ET<0]=0
    
    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['Prospero_ET'], mode='lines', name='Prospero_ET'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl2['Pristley_Taylor_ET'], mode='lines', name='Pristley_Taylor_ET'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl3['PenmanMontheithFAO_ET'], mode='lines', name='PenmanMontheithFAO_ET'))

   #fig.update_xaxes(rangeslider_visible=True)
    
    
    fig.update_layout(
        title='Compare EvapoTraspiration fluxes',
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        #xaxis_title="Date",
        yaxis_title="[$W\cdot m^{−2}$]",
        #legend_title="Date",
        font=dict(size=14))
    fig.show()


def compare_sim_obs(a,b):
    warnings.filterwarnings('ignore')
    ############################ GRAFICO 1 ########################################
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1)
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'GEOSPACE'
    kl.GEOSPACE[kl.GEOSPACE<0]=0
   
    
    kl2 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'SpikeII'
    kl2.SpikeII[kl2.SpikeII<0]=0
    
    
    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['GEOSPACE'], mode='lines', name='GEOSPACE'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl2['SpikeII'], mode='lines', name='SpikeII'))

   #fig.update_xaxes(rangeslider_visible=True)
    
    
    fig.update_layout(
        title='Compare GEOSPACE and SpikeII data',
        #xaxis_title="Date"
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        yaxis_title="$Evapotranspiration -[mm h^{-1}]$",
        #legend_title="Date",
        font=dict(size=16))
    fig.show()
    
    
def compare3(a,b,c):
    warnings.filterwarnings('ignore')
    ############################ GRAFICO 1 ########################################
    kl = pd.read_csv(a,skiprows=6,parse_dates=[1])
    kl = kl.drop(['Format'],axis=1)
    kl.columns.values[0] = 'Date'
    kl.columns.values[1] = 'GEOSPACE'
    kl.GEOSPACE[kl.GEOSPACE<0]=0
    
    kl2 = pd.read_csv(b,skiprows=6,parse_dates=[1])
    kl2 = kl2.drop(['Format'],axis=1)
    kl2.columns.values[0] = 'Date'
    kl2.columns.values[1] = 'Dataset_1'
    kl2.Dataset_1[kl2.Dataset_1<0]=0
    
    kl3 = pd.read_csv(c,skiprows=6,parse_dates=[1])
    kl3 = kl3.drop(['Format'],axis=1)
    kl3.columns.values[0] = 'Date'
    kl3.columns.values[1] = 'Dataset_2'
    kl3.Dataset_2[kl3.Dataset_2<0]=0
    
    
    fig = px.line()
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl['GEO-SPACE'], mode='lines', name='GEO-SPACE'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl2['Dataset_1'], mode='lines', name='Dataset_1'))
    fig.add_trace(go.Scatter(x=kl['Date'], y=kl3['Dataset_2'], mode='lines', name='Dataset_2'))
    
    fig.update_layout(
        title='Compare GEO-SPACE and two observation dataset',
        #xaxis_title="Date"
        font_family="Times New Roman",
        font_color="Black",
        title_font_family="Times New Roman",
        title_font_color="Black",
        #yaxis_title="$Evapotranspiration -[mm h^{-1}]$",
        #legend_title="Date",
        font=dict(size=18))
    fig.show()
    
    
    
    
def watercontent(ds,depth1,depth2,start_date,end_date,freq,lab):
    '''
    This function extracts the average water content considering all control volumes between the given depths.     
    
    - ds: dataframe in which the GEOSPACE output file is read. 
    
    - depth1: first depth starting from the surface, to be entered with the positive sign.
    - depth2: second depth starting from the surface, to be entered with the positive sign.
    
    - start_date: start date of the timeseries
                str,'dd-mm-yyyy hh:mm'
    - end_date: end date of the timeseries
                str,'dd-mm-yyyy hh:mm'
    More info (help(pd.date_range))
    
    - freq: frequency of the timeseries. 'H': hourly, 'D': daily, 'T': minutes
    - lab: label added at the end of the file name

        
    '''
    depth=ds.depth.values
    
    v1=ds.depth.sel(depth=[-depth1], method='nearest')
    V1=depth.tolist().index(v1)
    v2=ds.depth.sel(depth=[-depth2], method='nearest')
    V2=depth.tolist().index(v2)
    thetas=ds.theta.sel(depth=ds.depth.values[V2:V1], time=ds.time.values[:]).values
    theta=np.mean(thetas,axis=1)
    
    COLUMN_NAMES=['Datetime','vwc']
    df_hourly = pd.DataFrame(columns=COLUMN_NAMES)
    date_rng = pd.date_range(start=datainizio, end=datafine , freq=frequenza)
    df_hourly['Datetime']=date_rng
    df_hourly['vwc']=theta

    df_hourly.insert(loc=0, column='-', value=np.nan)
    ID_basin = 1
    write_timeseries_csv(df_hourly,'vwc_average_'+lab+'.csv',ID_basin)
