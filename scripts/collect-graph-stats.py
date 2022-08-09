#This file belongs to the "Seeking Critical Nodes in Digraphs\" project.
#The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs
#This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
#For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

import sys
import json
import argparse
from igraph import *
import traceback

METR=['betweenness-iterative', 'betweenness-standard', 'closeness_all-iterative', 'closeness_all-standard','cndp','degree_all-iterative','degree_all-standard','pagerank-iterative','pagerank-standard','bf','random']

parser=argparse.ArgumentParser(description='Collect information, metrics, analysis, run info about a given graph')
parser.add_argument('--name',help='graph name', required=True)
parser.add_argument('--filename', help='the graph filename; all the results of the analysis must be stored in the same directory, with same prefix name', required=True)
parser.add_argument('--graphtype', help='real/synthetic', choices=['real','sintetic'],required=True)
parser.add_argument('--graphfamily',choices=['roadnet','p2p','social','web','misc','ER','Lattice','RMAT','BA','WS'],required=True)
parser.add_argument('--fileformat',choices=['dimacs','rmat'],required=True)

args=parser.parse_args()
g=dict()
bf=args.filename
g['graph-name']=args.name
g['graph-type']=args.graphtype
g['graph-family']=args.graphfamily
print("[INFO]: reading file {}".format(bf))

try:
    bfj=bf.split(".out")[0]
    bfj='{}.json'.format(bfj)
    f=open(bfj,'r')
    j=json.load(f)

    g['root']=j['root']
    g['#nodes']=j['#nodes']
    g['#edges']=j['#edges']
    g['edges-list']=j['edges-list']
    g['global-properties']=j['global-properties']
    
except:
    
    try:
        with open(bf) as f:
            lines=f.readlines()
        root=0
        nodes=0
        edges=0
        edgeslist=[]
    
        if args.fileformat == 'dimacs':
            if 'p' not in lines[0]:
                print("Error: Invalid input-basefile {}".format(bf))
                print("DEBUG: file-header :{}".format(lines[0]))
                exit(1)
    
            a,nodes,edges,root=lines[0].replace('\n','').split(' ')
            root=root.replace('\n','')
    
            for i in range(1,len(lines)):
                a,src,dst=lines[i].replace('\n','').split(' ')
                src=int(src)
                dst=int(dst)
                edgeslist.append([src,dst])
        else:
            if '#Nodes:' not in lines[0]:
                print("Error: Invalid input-basefile {}".format(bf))
                print("DEBUG: file-header :{}".format(lines[0]))
                exit(1)
            a,nodes,b,edges=lines[0].replace('\n','').split(' ')  
            for i in range(1,len(lines)):
                src,dst=lines[i].replace('\n','').split(' ')
                src=int(src)
                dst=int(dst)
                edgeslist.append([src,dst])
            
        g['root']=int(root)
        g['#nodes']=int(nodes)
        g['#edges']=int(edges)
        g['edges-list']=edgeslist
    
        # Graph global properties
        print("[INFO]: computing global properties")
        gprop=dict()
        gr=Graph(n=int(nodes),edges=edgeslist,directed=True)
    
        print("[INFO]: computing reciprocity")
        gprop['reciprocity']=gr.reciprocity()
        print(gprop)
	
        g['global-properties']=gprop
    
    except Exception as e:
        print("Exception {}".format(e))
        print("Error: invalid input-basefile {}".format(bf))
        print(traceback.format_exc())
        exit(1)
    
bf=bf.split(".out")[0]
# METRICS
m_valid=True
for m in METR:
    m_valid=True
    print("[INFO]: reading {} results".format(m))
    try:
        with open("{}.{}.{}".format(bf,m,'out')) as f:
            lines=f.readlines()

    except Exception as e:
        print("Exception: {}".format(e))
        print("Skipping metric {}".format(m))
        m_valid=False
        continue
    data=dict()
    data['metric-name']=m
    
    data['computation-time']=0.0
    if 'standard' in m:
        time=float(lines[0].split(',')[0])
        data['computation-time']=time

    remlist=dict()
    for i in range(1,len(lines)):
        lines[i]=lines[i].replace('\n','')
        try:
            rand_chos=0
            if m =='bf':
                nodeid,conn=lines[i].split(',')
                time=0
            elif m == 'cndp':
                arr=lines[i].split(',')
                time=arr[0]
                nodeid=arr[1]
                rand_chos=arr[2]
                conn=arr[3]
            else:
                arr=lines[i].split(',')
                time=arr[0]
                nodeid=arr[1]
                conn=arr[2]
            remlist[i]={'k-size': i, 'nodeid' : int(nodeid), 'randomly_chosen': int(rand_chos), 'connectivity' : int(conn), 'time': float(time)}
        except:
            print("Metric {}: line {} malformed - skipping the rest of the file".format(m,i))
            print(lines[i])
            print(lines[i].split(',',6))
            i=len(lines)
            break
    data['removed-nodes']=remlist
    g[m]=data
if m_valid==True:
    if m == 'bf':
        initial_conn=int(lines[0].split(',')[1].replace('\n',''))
    else:
        initial_conn=int(lines[0].split(',')[2].replace('\n',''))
    g['start-connectivity']=initial_conn
 
try:
    #JSON WRITE
    bf='{}.{}'.format(bf,'json')
    with open(bf,'w+') as fout:
        json.dump(g,fout,indent=True)
except:
    print('Error: unable to dump to json {}'.format(bf))
    exit(1)





