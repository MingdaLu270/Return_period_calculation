
# coding: utf-8

# In[1]:


"""
This program downloads peak flow data from USGS Surface Data Portal for a USER_INPUT station and calculates the flow 
corresponding to (different) return period

This code is written in Python 3 format

Revision No: 05
Last Revised : 2020-02-16
"""

## Import the required Modules/Packages for obtaining the data from portal
import urllib.parse
import urllib.request

## Import the required Modules/Packages for calculating return period flow using Gamma Inverse Function
import math
import scipy.stats
import numpy as np
#import csv
from scipy.stats import gamma
from scipy.stats import invgamma


# In[2]:


## CELL-02
## Define a function for obtaining the peak flow data from USGS Surface Data Portal
def GetAnnualPeakFlowData_f(station_number,FolderName):
    """
    Input: Station Number, Folder Name
    Output: Peak Flow Values, Station Name
    """
    ## Building URLs
    var1 = {'site_no': station_number}
    part1 = 'https://nwis.waterdata.usgs.gov/nwis/peak?'
    part2 = '&agency_cd=USGS&format=rdb'
    link = (part1 + urllib.parse.urlencode(var1) + part2)
    print("The USGS Link is: \t")
    print (link)
    
    ## Opening the link & retrieving data
    response = urllib.request.urlopen(link)
    html = response.read()
    
    ## Assigning the location & Storing the Original Data
    
    #DataStore=FolderName + station_number + ".txt"
    with open(FolderName+'Data_' + station_number + '_raw'  + '.txt', 'wb') as f1:
        f1.write(html)
    f1.close
    
    ## Converts html from bytes class to str class
    html = html.decode()
    ## Splits the string by \n and converts list
    html2 = html.split('\r\n')
    
    ## To get the station name 
    line_no=0
    for line_no in range(len(html2)):
        ## Check if first six (use 0:7) characters is "#  USGS",
        if html2[line_no][0:7]=="#  USGS":
            station_name=html2[line_no][3:]
            break
        line_no+=1
    
    ## Define an empty string
    reqd_data = ''
    reqd_flow_list=[]
    for line in html2[74:]:
        ## Splits each line to col by tab separator
        cols = line.split('\t')
        if len(cols) == 1:
            continue
        ## Joins only date and peakflow
        ## cols[2] corresponds to Date of peak streamflow (format YYYY-MM-DD)
        ## cols[4] corresponds to Annual peak streamflow value in cfs
        newline = ','.join([cols[2],cols[4]])
        reqd_data += newline + '\n'
        reqd_flow_list.append((cols[4]))

    
    ## Converts reqd_data from str class to bytes class
    reqd_data = reqd_data.encode() 
    ## Saves the date and peakflow into a new file
    with open(FolderName+'Data_' + station_number + '_reqd'  + '.txt', 'wb') as f2:
        f2.write(reqd_data)
    f2.close
    print ('\n')
    print("Raw Data and Processed Data is stored in Results Folder.")
    
    ## Returns the peak flow data as list for calculation of return period
    return (reqd_flow_list,station_name)


# In[4]:


## CELL-03
## Main Code

station_number=input("Enter USGS Station Number of the Required Station (USGS Station Number/site_no) \t")
print('\t')
FolderName="./Results/"
peakflow_list_wb,station_name=GetAnnualPeakFlowData_f(station_number,FolderName)
print("\nThe station name is:", station_name,"\n")


# In[15]:


## CELL-04

## Enter the four years for carrying out the analysis
## Input data & analysis years

data_start_year=int(input("Enter the starting year of DATA PERIOD (excluding initial break period):"))
print('\t')
data_end_year=int(input("Enter the ending year of DATA PERIOD:"))
print('\t')
analysis_start_year=int(input("Enter the starting year of ANALYSIS PERIOD:"))
print('\t')
analysis_end_year=int(input("Enter the ending year of ANALYSIS PERIOD:"))
print('\t')


# In[20]:


## CELL-05

## WRITE YOUR CODE HERE
mw=30
## Time Step of Moving Window
ts_define=[1,5]
## Different Return Periods (in years)
rp_define=[100,500]

year_length=int(data_end_year-data_start_year+1)

peakflow_list=peakflow_list_wb[-year_length:]
data_length=len(peakflow_list)

length_rp=len(rp_define)

if data_length!=year_length:
    print('Error in length of data and chek weather it is continuous')
    
else:
    
    period_start=analysis_start_year-data_start_year
    period_end=(analysis_end_year-(data_start_year-1))-mw
    analysis_dur=period_end-period_start
    ## For looping for different time step
    for time_step in ts_define:
        row=1
        flow_matrix=np.zeros([2+math.floor(analysis_dur/time_step),len(rp_define)+2])
        
        ## To make the first row to hold the return period years
        
        for i in range(length_rp):
            flow_matrix[0][i+2]=int(rp_define[i])
            
        for j in range(period_start,period_end+1,time_step):
            sum=0
            sum_log=0
            for i in range(mw):
                sum_log=sum_log+math.log(int(peakflow_list[j+i]))
            mean=sum/mw
            mean_log=sum_log/mw
            diff_sum2=0
            diff_sum3=0

            ## Calculation of parameters
            for k in range(mw):
                diff_log=(math.log(int(peakflow_list[j+k]))-mean_log)
                diff_sum2=diff_sum2+diff_log**2
                diff_sum3=diff_sum3+diff_log**3
            stddev=math.sqrt(diff_sum2/(mw-1))
            skewness=mw/((mw-1)*(mw-2))*diff_sum3/stddev**3
            shape=4/skewness**2
            scale=stddev*skewness/2
            location=mean_log-shape*scale
            c=np.empty([period_end+1,8])
            col=2

            ##First two columns reserved from start year and end year
            flow_matrix[row][0]=analysis_start_year+(row-1)*time_step
            flow_matrix[row][1]=data_start_year-1+mw+j

            ## Calculation of return period flow
            m=0
            for rp in rp_define:
                Prob=(1-1/rp)
                if scale>0:
                    b=(gamma.ppf(Prob,shape,0,scale)-shape*scale)/skewness
                else:
                    b=(-gamma.ppf(1-Prob,shape,0,-scale)-shape*scale)/skewness
                LN_q=mean_log+b*skewness
                peak_flow=round(math.exp(LN_q))
                flow_matrix[row][col]=peak_flow
                file_path_name=FolderName+'RP_Flow_'+(station_number)+'_'+str(analysis_start_year)+'_'+str(analysis_end_year)+'_mw'+str(mw)+'_ts'+'{:02d}'.format(time_step)+'.txt'
                col=col+1
            row=row+1
        np.savetxt(file_path_name,flow_matrix,delimiter=',',fmt='%1.0f')
    print('Moving Average Method is completed for: ',station_name,'\n')
    print('REMEMBER: Peak flow calculated in cfs (cubic feet per second)')

