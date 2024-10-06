import numpy as np
import scipy as sci
import matplotlib as mlt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import scipy.integrate as s_i
import pandas as pd
con='n'
#Constants
G = 6.67e-11 #Nm^2/kg^2
m_sun = 2e30 #Kg
r_ac = 5.326e12 #m Distance b/w Alpha centauri stars
v_e = 30000 #ms^-1   <V_earth> around Sun
t_se =  3.156e7#Sec
link = 'https://github.com/Astrobuuz/2-Body-Simulation/tree/main'
#Scaling Constants
K1 = G*t_se*m_sun/(r_ac**2*v_e)
K2 = v_e*t_se/r_ac
#Defingng function for 2B motion
def two_b_eqn(w,t,G,m1,m2):
    r1 = w[:3]
    r2 = w[3:6]
    v1 = w[6:9]
    v2 = w[9:12]
    r = sci.linalg.norm(r2-r1) #Magnitude of distance between
    a1 = K1*m2*(r2-r1)/(r**3) #a1 = dv1/dt
    a2 = K1*m1*(r1-r2)/(r**3) #a2 = dv2/dt
    v1,v2 = K2*v1,K2*v2
    r_derivs = np.concatenate((v1,v2)) #dr1/dt,dr2/dt,dv1/dt,dv2/dt
    derivs = np.concatenate((r_derivs,a1,a2))
    return derivs
#Defining COM property
def com_pos(M1,R1,M2,R2):
    com = (M1*R1+M2*R2)/(M1+M2)
    return com

def com_vel(M1,V1,M2,V2):
    com_v = (M1*V1+M2*V2)/(M1+M2)
    return com_v
try:
    parameters = pd.read_csv('Saved_Parameters.csv')
except FileNotFoundError:
    print(f'File not found, download default file from github {link}  quiting program')
    quit()
except Exception:
    print(f'Something went wrong, Noclue what')
    quit()
choice=10
print("NOTE: All the data are in terms of Solar Mass, Orbital Velocity of Earth, Distance b/w Alpha Centauri, Time period of Earth Sun system[1 year], and should be entered likewise only for good performance.")
print(f"Constants and scaling factors:\nMass of Sun: {m_sun}\nOrbital Velocity of Earth: {v_e}\nDistance: {r_ac},\nTime Period: {t_se}\nForce Scaler: {K1}\nVelocity Scaler: {K2}")
print("===================================================================")
print("          Welcome to the 2 body Simulatior     ")
while choice != 5:
    print("===============================================================")
    print("1. Run Simulation")
    print("2. Export CSV solution")
    print("3. Import CSV and plot data")
    print("4. Edit the CSV")
    print("5. Exit")
    print("================================================================")
    if con=='y':
        choice=4
    else:
        choice = eval(input('Enter your choice (1/2/3/4/5): ')) 
    if choice==1:
        use_saved=input("Do you want to use saved parameters from csv(y/n)?:")
        if use_saved.lower()=='y':
            index_to_use = eval(input(f"{parameters.System_Name}\nEnter the index which you want to use:"))
            if index_to_use not in parameters.index:
                print("Index not found Using Alpha Centauri System")
                index_to_use = 0
            print(f"Using {parameters.System_Name[index_to_use]} system")
            m1 = parameters.m1[index_to_use]
            m2 = parameters.m2[index_to_use]
            r1 = np.array([parameters.r1_x[index_to_use],parameters.r1_y[index_to_use],parameters.r1_z[index_to_use]])
            r2 = np.array([parameters.r2_x[index_to_use],parameters.r2_y[index_to_use],parameters.r2_z[index_to_use]])
            v1 = np.array([parameters.v1_x[index_to_use],parameters.v1_y[index_to_use],parameters.v1_z[index_to_use]])
            v2 = np.array([parameters.v2_x[index_to_use],parameters.v2_y[index_to_use],parameters.v2_z[index_to_use]])
            com = com_pos(m1,r1,m2,r2)
            com_v = com_vel(m1, v1, m2, v2)
            fra = int(input('Enter number of datapoints[recommended 1500-2000]:'))
            time_span = np.linspace(0,float(parameters.timespan[index_to_use]),fra)
            body_1,body_2 = parameters.body_1[index_to_use],parameters.body_2[index_to_use]
        elif use_saved.lower()=='n':
            print("Using Custom Parameters:")
            body_1 = input('Enter the name of the first body:')
            body_2 = input('Enter name of the second body:')
            r1 = input(f"Enter the position of the {body_1}(x,y,z) separtated by commas:")
            r2 = input(f"Enter the position of the {body_2}(x,y,z) separtated by commas:")
            v1 = input(f"Enter the velocity of the {body_1}(vx,vy,vz) separtated by commas:")
            v2 = input(f"Enter the velocity of the {body_2}(vx,vy,vz) separtated by commas:")
            r1,r2 = '['+r1+']','['+r2+']'
            r1,r2 = eval(r1),eval(r2)
            v1,v2 = '['+v1+']','['+v2+']'
            v1,v2 = eval(v1),eval(v2)
            m1=float(input(f'Enter Mass of {body_1}:'))
            m2=float(input(f'Enter Mass {body_2}:'))
            timespan = int(input('Enter timespan:'))
            fra = int(input('Enter number of datapoints[recommended 1500-2000]:'))
            time_span = np.linspace(0,timespan,fra)
            save_to_csv = input('Would you like to csv (y/n):')
            if save_to_csv.lower() =='y':
                print(f'Saving to CSV')
                sys_name = input('Enter the name of your system:')
                parameters.loc[len(parameters)]=[sys_name,body_1,body_2,m1,m2,r1[0],r1[1],r1[2],r2[0],r2[1],r2[2],v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],timespan]
                print(f"This data will be saved for future use {parameters.tail(1)}")
                parameters.to_csv('Saved_Parameters.csv')
        init_para=np.array([r1,r2,v1,v2]).flatten()
        two_b_sol=s_i.odeint(two_b_eqn,init_para,time_span,args=(G,m1,m2))
        r1_sol,r2_sol=two_b_sol[:,:3],two_b_sol[:,3:6]
        fig=plt.figure()
        ax=fig.add_subplot(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('2B Problem Simulation')
        h1=[ax.scatter(r1_sol[0,0],r1_sol[0,1],r1_sol[0,2],c='b',s=80,label=body_1)]
        h2=[ax.scatter(r2_sol[0,0],r2_sol[0,1],r2_sol[0,2],c='r',s=80,label=body_2)]
        plt.legend()
        def anim(i):          
            h1[0].remove()
            h2[0].remove()
            ax.plot(r1_sol[:i,0],r1_sol[:i,1],r1_sol[:i,2],color = 'b')
            ax.plot(r2_sol[:i,0],r2_sol[:i,1],r2_sol[:i,2],color = 'r')
            h1[0]=ax.scatter(r1_sol[i,0],r1_sol[i,1],r1_sol[i,2],color="b",marker="o",s=80)
            h2[0]=ax.scatter(r2_sol[i,0],r2_sol[i,1],r2_sol[i,2],color="r",marker="o",s=80)
        ani = animation.FuncAnimation(fig,anim,frames=fra,interval=1)
        save=input("Do you want to  save the animation created(y/n)?:")
        if save.lower()=='y':
            name=input('Enter File name without extension:')
            name=name+'.mp4'
            print(f"Animation will be saved in the current directory as {name}")
            ani.save(name,writer='ffmpeg',dpi=300,fps=30)
            print("Animation saved successfully! :D")
        plt.show()
    elif choice==2:
        use_saved = input("Do you want to use saved parameters from the CSV? (y/n): ")
        if use_saved.lower() == 'y':
            index_to_use = eval(input(f"{parameters.System_Name}\nEnter the index which you want to use:"))
            if index_to_use not in parameters.index:
                print("Index not found Using Alpha Centauri System")
                index_to_use = 0
            else:
                print(f"Using {parameters.System_Name[index_to_use]} system")
                m1 = parameters.m1[index_to_use]
                m2 = parameters.m2[index_to_use]
                r1 = np.array([parameters.r1_x[index_to_use],parameters.r1_y[index_to_use],parameters.r1_z[index_to_use]])
                r2 = np.array([parameters.r2_x[index_to_use],parameters.r2_y[index_to_use],parameters.r2_z[index_to_use]])
                v1 = np.array([parameters.v1_x[index_to_use],parameters.v1_y[index_to_use],parameters.v1_z[index_to_use]])
                v2 = np.array([parameters.v2_x[index_to_use],parameters.v2_y[index_to_use],parameters.v2_z[index_to_use]])
                com = com_pos(m1,r1,m2,r2)
                com_v = com_vel(m1, v1, m2, v2)
                fra = int(input('Enter number of datapoints[recommended 1500-2000]:'))
                time_span = np.linspace(0,parameters.timespan,fra)
        elif use_saved.lower()=='n':
            print("Using Custom Parameters:")
            body_1 = input('Enter the name of the first body:')
            body_2 = input('Enter name of the second body:')
            r1 = input(f"Enter the position of the {body_1}(x,y,z) separtated by commas:")
            r2 = input(f"Enter the position of the {body_2}(x,y,z) separtated by commas:")
            v1 = input(f"Enter the velocity of the {body_1}(vx,vy,vz) separtated by commas:")
            v2 = input(f"Enter the velocity of the {body_2}(vx,vy,vz) separtated by commas:")
            r1,r2 = '['+r1+']','['+r2+']'
            r1,r2 = eval(r1),eval(r2)
            v1,v2 = '['+v1+']','['+v2+']'
            v1,v2 = eval(v1),eval(v2)
            m1=float(input(f'Enter Mass of {body_1}:'))
            m2=float(input(f'Enter Mass {body_2}:'))
            timespan = int(input('Enter timespan:'))
            fra = int(input('Enter number of datapoints[recommended 1500-2000]:'))
            time_span = np.linspace(0,timespan,fra)
            init_para = np.array([r1,r2,v1,v2]).flatten()
            two_b_sol = s_i.odeint(two_b_eqn,init_para,time_span,args=(G,m1,m2))
            r1_sol = two_b_sol[:,:3]
            r2_sol = two_b_sol[:,3:6]
            v1_sol = two_b_sol[:,6:9]
            v2_sol = two_b_sol[:,9:12]
            df = pd.DataFrame({
            'Time':time_span,
            'r1_x': r1_sol[:, 0],
            'r1_y': r1_sol[:, 1],
            'r1_z': r1_sol[:, 2],
            'r2_x': r2_sol[:, 0],
            'r2_y': r2_sol[:, 1],
            'r2_z': r2_sol[:, 2],
            'v1_x': v1_sol[:, 0],
            'v1_y': v1_sol[:, 1],
            'v1_z': v1_sol[:, 2],
            'v2_x': v2_sol[:, 0],
            'v2_y': v2_sol[:, 1],
            'v2_z': v2_sol[:, 2]
            })
            path_and_name = input('Enter the file path of and name without extension:')
            path_and_name = path_and_name+'.csv'
            df.to_csv(path_and_name)
            print(f'File saved as {path_and_name}')
    elif choice==3:
        path_to_csv = input('Enter the Address of the csv:\n')
        dtf = pd.read_csv(path_to_csv)
        #Taking Magnitudes of vector
        dtf['r1_mag'] = (dtf.r1_x**2+dtf.r1_y**2+dtf.r1_y**2)**(1/2)
        dtf['r2_mag'] = (dtf.r2_x**2+dtf.r2_y**2+dtf.r2_y**2)**(1/2)
        dtf['v1_mag'] = (dtf.v1_x**2+dtf.v1_y**2+dtf.v1_y**2)**(1/2)
        dtf['v2_mag'] = (dtf.v2_x**2+dtf.v2_y**2+dtf.v2_y**2)**(1/2)
        plot_choice = int(input('Do you want to\n1.Plot all the data\n2.Plot the distance from COM\n3.Plot the Velocity\n'))
        plt.xlabel('Time')
        r1_mag,r2_mag,v1_mag,v2_mag,t = dtf.r1_mag,dtf.r2_mag,dtf.v1_mag,dtf.v2_mag,dtf.Time
        if plot_choice == 1:
            plt.plot(r1_mag,t,label = 'r1')
            plt.plot(r2_mag,t,label = 'r2')
            plt.plot(v1_mag,t,label = 'v1')
            plt.plot(v2_mag,t,label = 'v2')
        elif plot_choice == 2:
            plt.plot(r1_mag,t,label = 'r1')
            plt.plot(r2_mag,t,label = 'r2')
        elif plot_choice == 3:
            plt.plot(v1_mag,t,label = 'v1')
            plt.plot(v2_mag,t,label = 'v2')
        else:print('Invalid Choice')
        plt.legend()
        plt.show()
    elif choice==4:
        print(f"Would you like to\n1.Edit a column\n2.Delete a data?:")
        edit_choice = eval(input(''))
        while edit_choice not in (1,2):
            edit_choice = input('Pls enter 1 or 2 only')    
        if edit_choice == 1:
            print(f"Which row would you like to change:\n{parameters.System_Name}")
            row_to_edit = eval(input("Enter row index:"))
            while row_to_edit not in parameters.index:
                print("Invalid index")
                row_to_edit = eval(input("Please enter a valid row index:"))
            print(f"Which column do you want to change:\n{parameters.columns}")
            column_to_edit = input("Enter column index:")
            while column_to_edit not in parameters.columns:
                print("Invalid Index")
                column_to_edit = eval(input("Please enter a valid column:"))
            if column_to_edit in ['System_Name','body_1','body_2']:
                new_value = input(f"The old value if {parameters.at[row_to_edit,column_to_edit]}, the new value should be:")
            else:
                new_value = float(input(f"The old value if {parameters.at[row_to_edit,column_to_edit]}, the new value should be:"))
                parameters.at[row_to_edit,column_to_edit]=new_value
        else:
            print(f'Which row would you like to delete?\n{parameters.System_Name}')
            row_to_delete = eval(input("Enter the row index:"))
            while row_to_delete not in parameters.index:
                row_to_delete = eval(input("Please enter a valid row index:"))
            parameters=parameters.drop(row_to_delete,axis='rows')
        print(f'Updated csv is:\n{parameters}')
        print(f"Should the changes be made to the file?(y/n):")
        update_to_csv=input('')
        if update_to_csv.lower()=='y':
            parameters.to_csv('Saved_Para_Edited.csv',index=False)
        con = input('Do you want to continue editing?:(y/n)')
