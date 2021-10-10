#importing libraries
from math import comb
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#Constructing Graph 
nodes=100 #number of nodes 
G_=nx.barabasi_albert_graph(nodes,2)
#G=nx.Graph()
#elist=[(0,2),(1,2),(1,3),(2,3),(2,4),(3,4)]
#G.add_edges_from(elist)

#define functions for objective functions for curve fitting
def objective1(x, a ,f):
	return a*x + f

def objective2(x, a ,b,f):
	return a*(x**b) + f


#taking data
num_walkers=50  #number of walkers
time=1000 #number of steps taken
walkers_pos_trac=np.zeros((time,num_walkers)) 
walkers_init_pos=np.random.randint(0,nodes,num_walkers) #initial position of a walkers
walkers_pos_trac[0,::]=walkers_init_pos



#defining a function that gives the next position of a walker 
def next_pos(pos):
    neig=G_.adj[pos]
    neig_list=[str(i) for i in neig.keys()]
    #print(neig_list)
    next_pos=np.random.choice(neig_list)
    #print(next_pos)
    return int(next_pos)

#a matrix containing positions of walkers with the evolution of time
#row represents position of a walker 
#column represents a particular walker
for j in range(num_walkers):
  for i in range(time):
    if i==time-1:
        break
    k=next_pos(walkers_pos_trac[i,j])
    walkers_pos_trac[i+1,j]=k
#nx.draw(G,with_labels=True)
#plt.show()

#numbers of walkers on a particular node at different time
#row represents the  descrete time
#Column represents the nodes
num_wal_nodes_=np.zeros((time,nodes))
for i in range(time):
    for j in range(0,nodes):
        m=walkers_pos_trac[i,::]
        c=np.count_nonzero(m==j)
        num_wal_nodes_[i,j]=c
#plotting a graph showing number of walkers on a particular node with time
plt.xlabel('Time(s)')
plt.ylabel('Number of walkers on a given node(5) ')
plt.plot(num_wal_nodes_[::,5])
plt.grid()
plt.title('Number of walkers on a given node at different time') 
plt.show()        
#calculating total number of walkers on a particular node from the whole simuation        
bin_dis=np.zeros((num_walkers+1,nodes))
for j in range(nodes):    
  for i in range(0,num_walkers+1):
    b=num_wal_nodes_[::,j]
    count=np.count_nonzero(b==i)
    bin_dis[i,j]=count
    
#calculating nodes having same degree 
dict_={v:d for v,d in G_.degree()}
set_=set(d for v,d in G_.degree())
uniq_deg=list(set_)
mat=-np.ones((max(set_)+1,nodes)) 
for j in range(len(uniq_deg)): 
  c=0
  deg=uniq_deg[j]
  for i in range(nodes):  
    k=dict_.get(i)
    if k==deg and k!=0:
      mat[deg,c]=i  
      c=c+1
      
#calculating average walkers on a particular nodes from whole simulation
avg_walk_node=np.zeros(max(set_)+1)
std_walk_node=np.zeros(max(set_)+1)
for i in range(max(set_)+1):
   mm=0
   cout=0
   std=0
   for j in range(nodes):
       if mat[i,j]==-1:
           break
       else:
           cout+=1
           aa=mat[i,j]
           mm+=np.mean(num_wal_nodes_[::,int(aa)])
           std+=np.std(num_wal_nodes_[::,int(aa)])
   if cout!=0:
     avg_walk_node[i]=mm/cout
     std_walk_node[i]=std/cout        
           
avg_walk_node_=[]
std_walk_node_=[]
for i,j in enumerate(avg_walk_node):
    if j!=0.0:
      avg_walk_node_.append([i,j]) 
for i,j in enumerate(std_walk_node):
    if j!=0.0:
      std_walk_node_.append([i,j])       

avg_walk_node_=np.array(avg_walk_node_)
std_walk_node_=np.array(std_walk_node_)


#plotting the graphs
x1=avg_walk_node_[::,0]
y1=avg_walk_node_[::,1]
popt, _ = curve_fit(objective1, x1, y1)
a,f = popt
y_line = objective1(x1, a, f)
plt.plot(x1,y_line,":",color='r',label='curve fitting($y=ax+c$)')
plt.scatter(x1,y1,label='Computed data points')
plt.grid()
plt.legend()
a=round(a,2)
c=round(f,2)
plt.text(max(set_)+2,1.6,'a='+str(a)+'\n c='+str(c))
plt.title('Curve fitting')
plt.xlabel("Degree of nodes")
plt.ylabel("Average number of walkers")
plt.show()

x2=std_walk_node_[::,0]
y2=std_walk_node_[::,1]
popt, _ = curve_fit(objective2, x2, y2)
a,b,f = popt
y_line = objective2(x2, a,b, f)
plt.plot(x2,y_line,":",color='r',label='curve fitting($y=ax^b+c$)')
plt.scatter(x2,y2,label='Computed data points')
plt.grid()
plt.legend()
a=round(a,2)
b=round(b,2)
c=round(f,2)
plt.text(max(set_)+2,1.6,'a='+str(a)+', b='+str(b)+'\n c='+str(c))
plt.title('Curve fitting')
plt.xlabel("Degree of nodes")
plt.ylabel("Standard deviation of walkers")
plt.show()
#occupation probability calculation of each node
#sum of all degrees of nodes
s=0
for v,d in G_.degree():
    s=s+d
occ_pb=np.zeros(nodes)
for i in range(nodes):
    occ_pb[i]=dict_[i]/float(s)  


#node with maximum degree
max_deg=max(set_)
node_max_deg=int(mat[max_deg,0])    

#histogram plotting of numerical simulation
bin_dis_=np.round_(bin_dis/time,3)
walkers_array=np.array([int(i) for i in range(0,num_walkers+1)])    
plt.bar(walkers_array,bin_dis_[::,node_max_deg],width=0.2)
plt.xlabel("ith walker \n on a given node "+str(node_max_deg))
plt.ylabel("Fraction of times \n ith walker is on a given node")
plt.grid()
plt.title("Binomial distribution from numerical simulation \n for node of maximum degree")
plt.show()    
    
#histogram plotting for theoretical calculation
z=np.zeros(num_walkers+1)
p=occ_pb[node_max_deg]
for i in range(num_walkers+1):
    z[i]=comb(50,i)*((p)**(i))*((1-p)**(50-i))
z=np.round_(z,3)
plt.bar(walkers_array,z,width=0.2) 
plt.grid()
plt.xlabel("ith walker on a \n given node "+str(node_max_deg))
plt.ylabel("Probability")
plt.title("Theoretically computed binomial distribution \nfor node of maximum degree")
plt.show()  
