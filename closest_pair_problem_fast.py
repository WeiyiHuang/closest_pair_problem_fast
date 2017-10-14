import numpy as np
import math

#select the points in specific section of dimension x
def locate_dimension_x(xl,xh,x,y):
    location_l=xl<x
    x=x[location_l]
    y=y[location_l]
    location_h=x<xh
    x=x[location_h]
    y=y[location_h]
    return x,y

#select the points in specific section of dimension u
def locate_dimension_y(yl,yh,x,y):
    location_l=yl<y
    x=x[location_l]
    y=y[location_l]
    location_h=y<yh
    x=x[location_h]
    y=y[location_h]
    return x,y

#select the points in specific section of dimension x and y, return the index of the points
def locate_xy(xl,xh,yl,yh,x,y):
    location_xl=x>xl
    location_xh=x<xh
    location_yl=y>yl
    location_yh=y<yh
    location=np.zeros((1,x.size),dtype=int)
    location=location.astype(bool)
    for i in range(0,location.size):
        location[0,i]=location_xh[i] and location_xl[i] and location_yh[i] and location_yl[i]
    p=np.argwhere(location[0]==True)
    if p.size==0:
        return -1
    else:
        return p[0][0]

#calculate Euclidean distance
def calculate_distance(x_1,x_2,y_1,y_2):
    return np.sqrt(np.square(x_1-x_2)+np.square(y_1-y_2))

#Divide and Conquer
'''divide points into the left half and the right half, perform recursive calls respectively on the two sides
to calculate the closest pair, compare them with the closest pair in the middle area '''
def divide_and_conquer(x,y,area_l,area_h):
    if x.size==1:
        closest_distance,x_a,x_b,y_a,y_b=bound_y_h-bound_y_l,-1,-1,-1,-1
    elif x.size==2:
        closest_distance=calculate_distance(x[0],x[1],y[0],y[1])
        x_a,x_b,y_a,y_b=x[0],x[1],y[0],y[1]
    else:
        area_mid=0.5*(area_h+area_l)
        l_x,l_y=locate_dimension_x(area_l,area_mid,x,y)
        r_x,r_y=locate_dimension_x(area_mid,area_h,x,y)
        if l_x.size>0 and r_x.size>0:
            distance_l,x_la,x_lb,y_la,y_lb=divide_and_conquer(l_x,l_y,area_l,area_mid)
            distance_r,x_ra,x_rb,y_ra,y_rb=divide_and_conquer(r_x,r_y,area_mid,area_h)
            temp=min(distance_l,distance_r)
            middle_x,middle_y=locate_dimension_x(area_mid-temp,area_mid+temp,x,y)
            distance_m,x_ma,x_mb,y_ma,y_mb=middle_area_calulation(middle_x,middle_y,area_mid,temp)
            if distance_m<distance_l and distance_m<distance_r:
                closest_distance=distance_m
                x_a,x_b,y_a,y_b=x_ma,x_mb,y_ma,y_mb
            else:
                if distance_l<distance_r:
                    closest_distance=distance_l
                    x_a,x_b,y_a,y_b=x_la,x_lb,y_la,y_lb
                else:
                    closest_distance=distance_r
                    x_a,x_b,y_a,y_b=x_ra,x_rb,y_ra,y_rb
        elif l_x.size==0:
            closest_distance,x_a,x_b,y_a,y_b=divide_and_conquer(r_x,r_y,area_mid,area_h)
        else:
            closest_distance,x_a,x_b,y_a,y_b=divide_and_conquer(l_x,l_y,area_l,area_mid)
    return closest_distance,x_a,x_b,y_a,y_b

#draw grid and select the closest pair in the middle area
def middle_area_calulation(x,y,area_mid,delta):
    grid=np.zeros(((int)(math.ceil((bound_y_h-bound_y_l)*2/delta)),4),dtype=int)
    for i in range(0,grid.shape[0]):
        for j in range(0,grid.shape[1]):
            grid[i,j]=locate_xy(area_mid-delta+j*delta/2,area_mid-delta+(j+1)*delta/2,bound_y_h-(i+1)*delta/2,bound_y_h-i*delta/2,x,y)
    d, x_ma,x_mb,y_ma,y_mb=bound_y_h-bound_y_l,-1,-1,-1,-1
    #for points in the first column, calculate the distance between itself and the 3 points in the third column
    j=0
    for i in range(0,grid.shape[0]):
        if grid[i,j]<>-1:
            m=2
            for n in range(-2,3):
                if ((i+n)>=0 and (i+n)<grid.shape[0] and grid[i+n,j+m]<>-1):
                    t=calculate_distance(x[grid[i,j]],x[grid[i+n,j+m]],y[grid[i,j]],y[grid[i+n,j+m]])
                    if t<d:
                        d=t
                        x_ma,x_mb,y_ma,y_mb=x[grid[i,j]],x[grid[i+n,j+m]],y[grid[i,j]],y[grid[i+n,j+m]]
    #for points in the second column, calculate the distance between itself and the 6 points in the last 2 column
    j=1
    for i in range(0,grid.shape[0]):
        if grid[i,j]<>-1:
            for m in range(1,3):
                for n in range(-2,3):
                    if ((i+n)>=0 and (i+n)<grid.shape[0] and grid[i+n,j+m]<>-1):
                        t=calculate_distance(x[grid[i,j]],x[grid[i+n,j+m]],y[grid[i,j]],y[grid[i+n,j+m]])
                        if t<d:
                            d=t
                            x_ma,x_mb,y_ma,y_mb=x[grid[i,j]],x[grid[i+n,j+m]],y[grid[i,j]],y[grid[i+n,j+m]]
    return d,x_ma,x_mb,y_ma,y_mb


#num of points
n=10
#define the initial section of the plane
bound_x_l,bound_y_l=0,0
bound_x_h,bound_y_h=1,1
#generate xs and ys randomly
x_0=np.random.rand(n)
y_0=np.random.rand(n)
print x_0
print y_0
#call the closest pari function
closest_distance,x_a,x_b,y_a,y_b=divide_and_conquer(x_0,y_0,bound_x_l,bound_x_h)
print closest_distance,x_a,x_b,y_a,y_b