#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 10:28:15 2021

CGTM Plotting Routines

@author: liam
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import all_possible_states as aps
import choose_path as chp
import matplotlib.colors as mcol
import CGTM_Calculations as cgc

class MidpointNormalize(mcol.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		mcol.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def plot_state_dist(pi_eq_file,ax):
    pi_eq = pd.read_csv("%s.csv"%pi_eq_file,index_col=0).T
    ax.bar(np.arange(0,len(aps.all_possible_states())),pi_eq.values[0])
    ax.set_ylabel("Prob of State")
    ax.set_xlabel("State")
    
    # plt.savefig("%sPi_Eq_%s.pdf"%(leaflet_in,kind))
    # plt.close()
    
def diff_plot(pi_eq_file,pi_raw_file, ax):
    pi = pd.read_csv("%s.csv"%pi_eq_file,index_col=0).T
    raw = pd.read_csv("%s.csv"%pi_raw_file,index_col=0).T
    dpi = pi.values[0] - raw.values[0]
    line = ax.bar(np.arange(0,len(aps.all_possible_states())), dpi)
    ax.set_ylabel(r"$\Delta$ Prob of State")
    ax.set_xlabel("State")
    
def plot_cgtm(TM_norm,ax):
    #TM_norm = pd.read_csv("%s.csv"%TM_norm_path,index_col=0,header=0)
    ax.pcolormesh(TM_norm,cmap="gist_earth_r")
    # ax.set_xlabel("State")
    # ax.set_ylabel("State")
    # plt.colorbar(label="Prob of Transition")
    # plt.savefig("%s_TMCG_%s.pdf"%(leaflet_in,kind))
    # plt.close()
    
def plot_state_traj(state,ax):
    states = pd.read_csv("%s.csv"%(state),index_col=0)
    
    if np.shape(states)[1]>1:
        s = ax.pcolormesh(states.T,cmap="ocean_r",vmin=0,vmax=230)
        ax.set_xlabel("Simulation")
        ax.set_ylabel("Frame")
        return s
    else:
        s = ax.pcolormesh(states,cmap="ocean_r",vmin=0,vmax=230)
        ax.set_ylabel("Frame")
        #ax.set_yticks(np.linspace(0,len(states)*.4,100))
        ax.set_xticks([])
        return s
    # plt.colorbar(label="State")
    # plt.savefig("%s_State_Change_Sims%s.pdf"%(leaflet_in,kind))
    # plt.close()

def plot_convergence(fl,ax):
    conv = pd.read_csv(fl,index_col=0)#[::-1]
    ax.plot(conv.index*.4, conv,"o-")
    ax.plot([0, ((conv.index+1)*.4)[-1]],[np.mean(conv),np.mean(conv)],'r--',lw=2)
    ax.plot([0, ((conv.index+1)*.4)[-1]],[1,1],'k--',lw=2)
    # ax.set_yscale('log')
    #ax.set_ylim(.75,1.25)
    ax.set_xlabel(r"$\tau$(ns)")

def plot_sigConverge(sigSU, sigSL,kind):
    sim_list = []
    if len(sigSU) <= 10:
        sim_list = [1,2,3,4,5,6,7,8,9,10]
    else:
        sim_list = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]
    plt.plot(np.linspace(0,50,len(sigSU)),sigSU,"o--",label="Outer Leaflet")
    plt.plot(np.linspace(0,50,len(sigSU)),sigSL,"+--",label="Inner Leaflet")
    plt.ylabel(r"$\sigma_{\pi}$")
    plt.xlabel("Simulations Used")
    plt.legend()
    plt.show()
    # plt.savefig("Convergence_%s.pdf"%kind)
    # plt.close()

def Ternary_Heat_Map(leaflet_in,fl_name,leaflet_in2=None,ax=None,out=None,initial=None):

    import matplotlib.tri as tri

    def plot_ticks(start, stop, tick, n, offset=(.0, .0)):
        r = np.linspace(0, 1, n+1)
        x = start[0] * (1 - r) + stop[0] * r
        x = np.vstack((x, x + tick[0]))
        y = start[1] * (1 - r) + stop[1] * r
        y = np.vstack((y, y + tick[1]))
        ax.plot(x, y, 'k', lw=1)
        
        # add tick labels
        for xx, yy, rr in zip(x[1], y[1], r):
            ax.text(xx+offset[0], yy+offset[1], "{:.2}".format(rr))
            
            
    '''
    get raw and get hist, run raw, run pi eq could be inputs instead of 
    functions to collect data....
    '''        
    def get_raw(leaflet_in):
        states = pd.read_csv("%s/%s.csv"%(chp.choose_path()[1],leaflet_in),index_col=0).T   
        return states       
    def get_hist(states):
        hist,edge = np.histogram(states,bins=len(aps.all_possible_states()),normed=None,range=(0,len(aps.all_possible_states())))
        return hist/hist.sum(),edge
    
    def run_pi_eq(state):
        states_r = get_raw(state).values[0]
        states = np.asarray(aps.all_possible_states())
        return  states,states_r 
    
    def run_raw(state):
        ## Should be removed...
        states_r = get_raw(state)
        hist ,edge = get_hist(states_r)
        states = np.asarray(aps.all_possible_states())
        return states, hist
    
    def run_ternary(state,fl_name,ax,out,initial):
        n = 4
        tick_size = 0.1
        margin = 0.05
        norm2 = MidpointNormalize(0,0.06,0.03)

        # define corners of triangle    
        left = np.r_[0, 0]
        right = np.r_[1, 0]
        top = np.r_[0.5,  np.sqrt(2.3)*0.576]
        triangle = np.c_[left, right, top, left]
        
        # define vectors for ticks
        bottom_tick = 0.8264*tick_size * (right - top) / n
        right_tick = 0.8264*tick_size * (top - left) / n
        left_tick = 0.8264*tick_size * (left - right) / n
        
        hist = 0
        states = 0

        # elif pi_eq == True:
        state = pd.read_csv("%s/%s.csv"%(chp.choose_path()[1],state),index_col=0).T
        states, hist = aps.all_possible_states(), state  

        # Note that the ordering from start to stop is important for the tick labels
        plot_ticks(right, left, bottom_tick, n, offset=(0, -0.06))
        plot_ticks(left, top, left_tick, n, offset=(-0.12, -0.0))
        plot_ticks(top, right, right_tick, n,offset=(0,.01))
        
        a = states[:,0]
        b = states[:,1]
        c = states[:,2]
        
        # # values is stored in the last column
        v = hist.values[0]
        
        # # translate the data to cartesian corrds
        x = 0.5 * ( 2.*b+c ) / ( a+b+c )
        y = 0.5*np.sqrt(3) * c / (a+b+c)
        
        
        # # create a triangulation out of these points
        T = tri.Triangulation(x,y)
        
        # # plot the contour
        if out == None:
            ax.tricontourf(x,y,T.triangles,v,cmap='RdBu_r',norm=norm2,levels=100)
        else:
            out = ax.tricontourf(x,y,T.triangles,v,cmap='RdBu_r',norm=norm2,levels=100)
        
        # create the grid
        corners = np.array([[0, 0], [1, 0], [0.5,  np.sqrt(2.3)*0.576]])
        triangle = tri.Triangulation(corners[:, 0], corners[:, 1])
        
        # creating the grid
        refiner = tri.UniformTriRefiner(triangle)
        trimesh = refiner.refine_triangulation(subdiv=3)
        tern = ax.triplot(trimesh,'--',color='grey')
        
        if initial is not None:
            hist_mod = pd.read_csv("%s/%s.csv"%(chp.choose_path()[1],initial),index_col=0,header=0)
    
            a1,b1,c1 = hist_mod.iloc[:,0], hist_mod.iloc[:,1], hist_mod.iloc[:,2]
            a1 = a1[a1>0].values
            b1 = b1[b1>0].values
            c1 = c1[c1>0].values
            x = 0.5 * ( 2.*b1+c1 ) / ( a1+b1+c1 )
            y = 0.5*np.sqrt(3) * c1 / (a1+b1+c1)
            # points = np.array([a,b,c]).T
            tri_point = tri.Triangulation(x,y)
            refiner = tri.UniformTriRefiner(tri_point)
            trimesh = refiner.refine_triangulation(subdiv=2)
            tern = ax.triplot(tri_point,'*',color='black')

        ax.set_axis_off()
        if out == None:
            return
        else:
            sm = plt.cm.ScalarMappable(norm=norm2, cmap = out.cmap)
            return out,sm        
    def run_ternary_diff(statem,staten,ax,out):
        import matplotlib.colors as mcolors
        n = 4
        tick_size = 0.1
        margin = 0.05
        
        # define corners of triangle    
        left = np.r_[0, 0]
        right = np.r_[1, 0]
        top = np.r_[0.5,  np.sqrt(2.3)*0.576]
        triangle = np.c_[left, right, top, left]
        
        # define vectors for ticks
        bottom_tick = 0.8264*tick_size * (right - top) / n
        right_tick = 0.8264*tick_size * (top - left) / n
        left_tick = 0.8264*tick_size * (left - right) / n
        
        # state = "pi_sl"#"StatesL_I"
        hist1,hist2 = 0,0
        states1,states2 = 0,0
        # if raw == True:
        #     states, hist = run_raw(state)    
        
        # elif pi_eq == True:
        state1 = pd.read_csv("%s/%s.csv"%(chp.choose_path()[1],statem),index_col=0).T
        states, hist1 = aps.all_possible_states(), state1  

        state2 = pd.read_csv("%s/%s.csv"%(chp.choose_path()[1],staten),index_col=0).T
        states, hist2 = aps.all_possible_states(), state2  

        #Define twin axis
        # Note that the ordering from start to stop is important for the tick labels
        plot_ticks(right, left, bottom_tick, n, offset=(0, -0.06))
        plot_ticks(left, top, left_tick, n, offset=(-0.15, -0.0))
        plot_ticks(top, right, right_tick, n,offset=(0,.01))
        
        a = states[:,0]
        b = states[:,1]
        c = states[:,2]
        
        # # values is stored in the last column
        v = (hist1 - hist2) #/ np.sum((hist1.T['0']- hist2.T['0'])) 
        norm2 = MidpointNormalize(midpoint=0,vmin=-1E-1,vmax=1E-1)#(vmin=v.values[0].min(),vmax=v.values[0].max(),midpoint=(v.values[0].max()+v.values[0].min())/2)
        
        # # translate the data to cartesian corrds
        x = 0.5 * ( 2.*b+c ) / ( a+b+c )
        y = 0.5*np.sqrt(3) * c / (a+b+c)
        
        
        T = tri.Triangulation(x,y)

        out = ax.tricontourf(x,y,T.triangles,v.T['0'],cmap='PuOr',norm=norm2,extend='both',levels=100)
        
        # create the grid
        corners = np.array([[0, 0], [1, 0], [0.5,  np.sqrt(2.3)*0.576]])
        triangle = tri.Triangulation(corners[:, 0], corners[:, 1])
        
        # creating the grid
        refiner = tri.UniformTriRefiner(triangle)
        trimesh = refiner.refine_triangulation(subdiv=3)
        
        #plotting the mesh and caliberate the axis
        ax.triplot(trimesh,'k--')
        ax.set_axis_off()

        if out == None:
            return
        else:
            sm = plt.cm.ScalarMappable(norm=norm2, cmap = out.cmap)
            return out,sm

    def run_ternary_iterate(state,ax,out):
            import matplotlib.colors as mcolors
            n = 4
            tick_size = 0.1
            margin = 0.05
            norm2 = 0
            norm2 = MidpointNormalize(0,0.12,0.06)

    
            # define corners of triangle    
            left = np.r_[0, 0]
            right = np.r_[1, 0]
            top = np.r_[0.5,  np.sqrt(2.3)*0.576]
            triangle = np.c_[left, right, top, left]
            
            # define vectors for ticks
            bottom_tick = 0.8264*tick_size * (right - top) / n
            right_tick = 0.8264*tick_size * (top - left) / n
            left_tick = 0.8264*tick_size * (left - right) / n
            
            # state = "pi_sl"#"StatesL_I"
            hist = 0
            states = 0
            # if raw == True:
            #     states, hist = run_raw(state)    
            
            # elif pi_eq == True:
            # state = pd.read_csv("%s/%s.csv"%(chp.choose_path()[1],state),index_col=0).T
            states, hist = aps.all_possible_states(), state  
    
            # Note that the ordering from start to stop is important for the tick labels
            plot_ticks(right, left, bottom_tick, n, offset=(0, -0.06))
            plot_ticks(left, top, left_tick, n, offset=(-0.12, -0.0))
            plot_ticks(top, right, right_tick, n,offset=(0,.01))

            a = states[:,0]
            b = states[:,1]
            c = states[:,2]
            
            # # values is stored in the last column
            v = hist.values#[0]#["0"]
            
            # # translate the data to cartesian corrds
            x = 0.5 * ( 2.*b+c ) / ( a+b+c )
            y = 0.5*np.sqrt(3) * c / (a+b+c)
            
            
            # # create a triangulation out of these points
            T = tri.Triangulation(x,y)
            
            # # plot the contour

            out = ax.tricontourf(x,y,T.triangles,v,cmap='RdBu_r',norm=norm2,levels=100,extend="both")
            
            # create the grid
            corners = np.array([[0, 0], [1, 0], [0.5,  np.sqrt(2.3)*0.576]])
            triangle = tri.Triangulation(corners[:, 0], corners[:, 1])
            
            # creating the grid
            refiner = tri.UniformTriRefiner(triangle)
            trimesh = refiner.refine_triangulation(subdiv=3)
            tern = ax.triplot(trimesh,'--',color='grey')

            if out == None:
                return
            else:
                sm = plt.cm.ScalarMappable(norm=norm2, cmap = out.cmap)
                return out,sm       
    def run_ternary_diff_iter(state1,state2,ax,out):
            import matplotlib.colors as mcolors
            n = 4
            tick_size = 0.1
            margin = 0.05
            
            # define corners of triangle    
            left = np.r_[0, 0]
            right = np.r_[1, 0]
            top = np.r_[0.5,  np.sqrt(2.3)*0.576]
            triangle = np.c_[left, right, top, left]
            
            # define vectors for ticks
            bottom_tick = 0.8264*tick_size * (right - top) / n
            right_tick = 0.8264*tick_size * (top - left) / n
            left_tick = 0.8264*tick_size * (left - right) / n
            
            # state = "pi_sl"#"StatesL_I"
            hist1,hist2 = 0,0
            states1,states2 = 0,0
            # if raw == True:
            #     states, hist = run_raw(state)    
            
            # elif pi_eq == True:
            states, hist1 = aps.all_possible_states(), state1  
    
            states, hist2 = aps.all_possible_states(), state2  
    
            #Define twin axis
            # Note that the ordering from start to stop is important for the tick labels
            plot_ticks(right, left, bottom_tick, n, offset=(0, -0.06))
            plot_ticks(left, top, left_tick, n, offset=(-0.15, -0.0))
            plot_ticks(top, right, right_tick, n,offset=(0,.01))
            
            a = states[:,0]
            b = states[:,1]
            c = states[:,2]
            
            # # values is stored in the last column
            v = (hist1 - hist2) #/ np.sum((hist1.T['0']- hist2.T['0'])) 
            norm2 = MidpointNormalize(midpoint=0,vmin=-1E-1,vmax=1E-1)#(vmin=v.values[0].min(),vmax=v.values[0].max(),midpoint=(v.values[0].max()+v.values[0].min())/2)
            
            # # translate the data to cartesian corrds
            x = 0.5 * ( 2.*b+c ) / ( a+b+c )
            y = 0.5*np.sqrt(3) * c / (a+b+c)
            
            
            T = tri.Triangulation(x,y)
    
            out = ax.tricontourf(x,y,T.triangles,v.T['0'],cmap='PuOr',norm=norm2,extend='both',levels=100)
            
            # create the grid
            corners = np.array([[0, 0], [1, 0], [0.5,  np.sqrt(2.3)*0.576]])
            triangle = tri.Triangulation(corners[:, 0], corners[:, 1])
            
            # creating the grid
            refiner = tri.UniformTriRefiner(triangle)
            trimesh = refiner.refine_triangulation(subdiv=3)
            
            #plotting the mesh and caliberate the axis
            ax.triplot(trimesh,'k--')
            ax.set_axis_off()

    
            if out == None:
                return
            else:
                sm = plt.cm.ScalarMappable(norm=norm2, cmap = out.cmap)
                return out,sm
    if isinstance(leaflet_in, str) == True:
        if leaflet_in2 is not None :
            return run_ternary_diff(leaflet_in,leaflet_in2,ax,out)
        else:
            return run_ternary(leaflet_in,fl_name,ax,out,initial)
    elif isinstance(leaflet_in, pd.core.series.Series) or isinstance(leaflet_in, pd.core.frame.DataFrame) == True:
        if leaflet_in2 is None:
            return run_ternary_iterate(leaflet_in,ax,out)  
        else:
            return run_ternary_diff_iter(leaflet_in,leaflet_in2,ax,out)  

    
def network_plot(ax,leaflet=None, kind=None, act=None,inds=False):
    import all_possible_states as aps
    import networkx as nx
    import numpy as np
    from matplotlib import cm
    import itertools
    
    
    def rescale(l,newmin,newmax):
        arr = list(l)
        return list(np.array([(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]))
        
    
    
    list_o_list = [np.linspace(0,230,231),np.linspace(0,230,231)]
    states = list(itertools.product(*list_o_list))
    states = np.arange(0,len(aps.all_possible_states()))
    pos_states = aps.all_possible_states()
    Q = 0
    
    if act == None:
        Q = cgc.CGTM_Calculations(leaflet,1,kind,).build_CGTM()[-1]
    else:
        Q = cgc.CGTM_Calculations("",1,"cg",act,"short").build_CGTM()[-1]
    G = nx.MultiDiGraph()
    G2 = None
    # edge_labels={}
    #pos = []
    # labs = {}
    # k = 0
    for i, origin_state in enumerate(states):
        for j, destination_state in enumerate(states):
            rate = Q[i][j]
            if rate > 0:
                # k = k + 1
                #pos.append([i,j])
                G.add_edge(origin_state,
                    destination_state,
                    weight=rate,
                    label="%i->%i"%(i,j)
                )
                # labs[k] ="%i->%i"%(i,j)

                # edge_labels[(origin_state, destination_state)] = label="{:.02f}".format(rate)
    pos1 = {}
    for p in G.nodes(data=True):
        pos1[p[0]] = (pos_states[p[0]][2],pos_states[p[0]][1])
    nx.set_node_attributes(G,pos1,"coord")
    base_size = 1
    # labs2 = {}
    if inds is not None:
        G2 = nx.MultiDiGraph()
        Q2 = cgc.CGTM_Calculations("",1,"cg",act,"long").build_CGTM()[-1]

        # edge_labels={}
        #pos = []
        # k = 0
        for i, origin_state in enumerate(states):
            for j, destination_state in enumerate(states):
                rate = Q2[i][j]
                if rate > 0:
                    # k = k + 1
                    #pos.append([i,j])
                    G2.add_edge(origin_state,
                        destination_state,
                        weight=rate,
                        label="%i->%i"%(i,j)#label="{:.02f},{:.02f},{:.02f}".format(pos_states[int(j)][0],pos_states[int(j)][1],pos_states[int(j)][2])
                    )
                    # labs2[k]="%i->%i"%(i,j)
        pos2 = {}
        for p in G2.nodes(data=True):
            pos2[p[0]] = (pos_states[p[0]][2],pos_states[p[0]][1])
    
    pos = {**pos1,**pos2}
    G3 = nx.compose(G, G2)
    # pos3 = nx.drawing.nx_pydot.graphviz_layout(G3, prog="circo")
    #nx.draw(G3,pos=pos,with_labels=True,ax=ax)
    # nx.draw_networkx_nodes(G3,pos1,nodelist=G.nodes(),node_color='b',ax=ax,alpha=.25)
    # nx.draw_networkx_nodes(G3,pos2,nodelist=G2.nodes(),ax=ax,node_color='r',alpha=.25)
    # nx.draw_networkx_edges(G3,pos1,edgelist=G.edges(),ax=ax,edge_color="b",alpha=.25)
    # nx.draw_networkx_edges(G3,pos2,edgelist=G2.edges(),ax=ax,edge_color="r",alpha=.25)
    # nx.draw_networkx_labels(G3,pos,ax=ax)
    # nx.draw_networkx_labels(G3,pos3,labels=labs2,font_size=8,ax=ax)
    # nx.draw_networkx_labels(G3,pos3,labels=labs,font_size=8,ax=ax)

    ####
    #Old
    ####
    # graph_colormap = cm.get_cmap('Reds', 12)
    # # node color varies with Degree
    # c = rescale([G.degree(v) for v in G],0.0,0.9) 
    # c = [graph_colormap(i) for i in c]
    # # node size varies with betweeness centrality - map to range [10,100] 
    G4 = nx.DiGraph(G)
    eigen_centrality = nx.eigenvector_centrality(G4, max_iter=1000)
    bc = eigen_centrality.copy()#nx.betweenness_centrality(G) # betweeness centrality
    s1 =  rescale([v for v in bc.values()],500,2500)
    # edge width shows 1-weight to convert cost back to strength of interaction 
    ew1 = rescale([float(G[u][v][0]['weight']) for u,v,null in G.edges],0.1,4)
    # edge color also shows weight
    # ec = rescale([float(G[u][v][0]['weight']) for u,v,null in G.edges],0.1,1)
    # ec = [graph_colormap(i) for i in ec]
    G5 = nx.DiGraph(G2)
    eigen_centrality = nx.eigenvector_centrality(G5, max_iter=1000)
    bc = eigen_centrality.copy()#nx.betweenness_centrality(G) # betweeness centrality
    s2 =  rescale([v for v in bc.values()],500,2500)
    ew2 = rescale([float(G2[u][v][0]['weight']) for u,v,null in G2.edges],0.1,4)

    
    nx.draw_networkx_nodes(G3,pos1,nodelist=G.nodes(),node_color='b',node_size=s1,ax=ax,alpha=.25)
    nx.draw_networkx_nodes(G3,pos2,nodelist=G2.nodes(),ax=ax,node_color='r',node_size=s2,alpha=.25)
    nx.draw_networkx_edges(G3,pos1,edgelist=G.edges(),ax=ax,edge_color="b",width=ew1,alpha=.25)
    nx.draw_networkx_edges(G3,pos2,edgelist=G2.edges(),ax=ax,edge_color="r",width=ew2,alpha=.25)
    nx.draw_networkx_labels(G3,pos,ax=ax)
    
    # pos = nx.drawing.nx_pydot.graphviz_layout(G, prog="circo")
    # nx.draw_networkx(G, pos=pos, with_labels=True, node_color=c, node_size=s,edge_color= ec,width=ew,
    #               font_color='white',font_weight='bold',font_size='9',ax=ax)
