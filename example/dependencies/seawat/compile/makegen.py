import os, sys, glob, re, shutil, glob
import networkx as nx
import matplotlib.pyplot as plt
import copy
from collections import Counter

def getFortFiles(srcDir,srcExt,srcExcl):

    file_list = []
    for d in srcDir:
        for e in srcExt:
            f = os.path.join(d,"*"+e)
            file_list.extend(glob.glob(f))

    # remove
    tmp = copy.deepcopy(file_list)
    for f in tmp:
        for ef in srcExcl:
            if ef == f:
                print "Excluding %s..."%ef
                file_list.remove(f)
    # check
    tmp = []
    for f in file_list:
        tmp.append(os.path.basename(f))
    d = [k for k, v in Counter(tmp).iteritems() if v > 1]
    if len(d) > 0:
        for fd in d:
            for f in file_list:
                if fd == os.path.basename(f):
                    print "Duplicate: %s"%f 
        raise Exception("Duplicate files found:")

    return file_list  
    
def createMod2File(file_list):
    mod_dict = {}
    for f in file_list:
        with open(f,'r') as fr:
            s = fr.read()
        lst = s.lower().split('\n')
        for s in lst:  
            if re.search(r'\bend module\b', s) == None:
                if re.search(r'\bmodule\b', s) != None:
                    s = s.strip()
                    if s[0:6] == 'module':
                        mod_dict[s.split()[1]] = os.path.basename(f)
    return mod_dict
    
def createFile2Use(file_list):
    use_dict = {}
    for f in file_list:
        fb = os.path.basename(f)
        useList = []
        with open(f,'r') as fr:
            s = fr.read()
        lst = s.lower().split('\n')
        for s in lst:  
            if re.search(r'\buse\b', s) != None:
                s = s.strip()
                if s[0:3] == 'use':
                    ss = s.split()[1].strip()
                    i = ss.find(",")
                    if i > -1:
                       ss = ss[:i]
                    useList.append(ss)
        if useList != []: use_dict[fb] = list(set(useList))
    return use_dict 
              
def createDiGraph(mod_dict,use_dict):
    
    G=nx.DiGraph()
    for f in use_dict:
        for mod in use_dict[f]:
            try:  
                fm = mod_dict[mod]
                if f != fm:
                    #print "%s: %s --> %s"%(mod,mod_dict[mod],f) 
                    if not G.has_node(fm):
                        G.add_node(fm)
                    if not G.has_node(f):
                        G.add_node(f)
                    if not G.has_edge(fm,f):
                        #print "edge: %s --> %s"%(fm,f)
                        G.add_edge(fm,f)
            except:
                print "Warning: module %s not found!"%mod
    return G

def saveGraph(G,cyc,nodsize,f):
    for node in G.nodes():
        if node in sum(cyc,[]):    
            G.node[node]['category'] = 'cyclic'
        else:
            G.node[node]['category'] = 'non-cyclic'
    color_map = {'cyclic':'r', 'non-cyclic':'g'}           
    nx.draw(G,with_labels=True,arrows=True,node_size=nodsize,node_color=[color_map[G.node[node]['category']] for node in G]) 
    plt.savefig(f)

if __name__ == "__main__":

    #######################################################
    bindir   = "/home/jengelen/seawat/bin/"
    program  = "seawat_svn317"
    f90flags = "-O3 -fpp -DOSD_OS_LINUX -DIFORT -DPKSMPI -DCLIP -DDEFLATION -DVTK -heap-arrays 10000000 -assume buffered_io" 
    f90      = "mpiifort"
    srcInDir  = ["/home/jengelen/seawat/imod-wq_pks/source/",
                 "/home/jengelen/seawat/imod-wq_pks/source/deltares_pck/",
                 "/home/jengelen/seawat/imod-wq_pks/source/deltares_utl/",
                 "/home/jengelen/seawat/imod-wq_pks/source/vtkfortran"]
    srcExcl   = ["/home/jengelen/seawat/imod-wq_pks/source/mt_fmi5.for",
                 "/home/jengelen/seawat/imod-wq_pks/source/gmg1.f",
                 "/home/jengelen/seawat/imod-wq_pks/source/vtkfortran/write_vtu.f90",
                 "/home/jengelen/seawat/imod-wq_pks/source/utl7u1vtu2.f"]
    srcOutDir = "./src/"
    srcExt    = ['.f','.F','.for','.FOR','.f90','.F90']
    hdrExt    = ['.inc','.h','.fi','.gin','.com']
    overwrite = True

    # reals as doubles
    if False:
        program  = "seawatr8"
        f90flags = "-O3 -fpp -DOSD_OS_LINUX -DIFORT -DPKSMPI -DCLIP -DDEFLATION -DVTK -DPKSDOUBLE -heap-arrays 1000000000 -assume buffered_io -real-size 64 -align dcommons" 
        srcOutDir = "./srcr8/"
    # reals as doubles debug
    if False:
        program  = "seawatr8-deb"
        f90flags = "-O0 -g -traceback -check bounds -debug all -fpp -DOSD_OS_LINUX -DIFORT -DPKSMPI -DCLIP -DDEFLATION -DPKSDOUBLE -heap-arrays 10000000 -assume buffered_io -real-size 64 -align dcommons"
        srcOutDir = "./srcr8_deb/"
        srcInDir = srcInDir[:-1]
    # debug
    if False:
        program  = "seawat_svn317-deb"
        f90flags = "-O0 -g -traceback -check bounds -debug all -fpp -DOSD_OS_LINUX -DIFORT -DPKSMPI -DCLIP -heap-arrays 10000000 -assume buffered_io" 
        srcOutDir = "./src_debug/"
        srcInDir = srcInDir[:-1]

    # openmp
    if False:
        program  = "seawat-omp"
        f90flags = "-O3 -fpp -openmp -DOSD_OS_LINUX -DIFORT -DPKSMPI -DCLIP -DPKSOMP -heap-arrays 10000000 -assume buffered_io" 
        srcOutDir = "./src_omp/"

    # scalasca
    scorep = False
    if False:
        scorep = True 
        f90      = "scorep mpiifort"
        program  = "seawat-sca"
        f90flags = "-O3 -fpp -DOSD_OS_LINUX -DIFORT -DPKSMPI -DCLIP -heap-arrays 10000000 -assume buffered_io" 
        srcOutDir = "./src_sca/"

    gperf = False
    if gperf:
        f90flags = f90flags+" -g"
        program  = "seawat-gperf"
        srcOutDir = "./src_gperf/"


    #######################################################

    # get entire list of all source files
    file_list = getFortFiles(srcInDir,srcExt+hdrExt,srcExcl)

    # create dictionary with modules and associated files
    mod_dict = createMod2File(file_list)       
    use_dict = createFile2Use(file_list)

    # create output directory
    if not os.path.isdir(srcOutDir): os.mkdir(srcOutDir)

    # copy the files
    print "Copying files..."
    i = 0
    for f in file_list:
        if overwrite:
            shutil.copy(f, srcOutDir) 
        else:
            # check is file already esists
            ft = os.path.join(srcOutDir,os.path.basename(f))
            if not os.path.isfile(ft):
                shutil.copy(f, srcOutDir) 
                print "%s --> %s"%(os.path.basename(f),srcOutDir) 

        #print "%s --> %s"%(os.path.basename(f),srcOutDir)
        file_list[i] = os.path.basename(f); i+= 1
        
    # rename .for file is case of scalasca
    if scorep:
        tmp1 = glob.glob(os.path.join(srcOutDir,'*.for'))
        tmp2 = glob.glob(os.path.join(srcOutDir,'*.FOR'))
        for f in tmp1+tmp2:
            fn = os.path.splitext(f)[0]+'.f'
            os.rename(f,fn)

    # create the graph
    G = createDiGraph(mod_dict,use_dict)   
 
    # check for cyclic dependencies
    cyc = list(nx.simple_cycles(G))
 
    # draw the graph
    # saveGraph(G,cyc,3000,"graph.pdf")
 
    # raise for cyclic dependency
    if cyc != []:
        raise Exception("Error, cyclic dependency found!")
    
    file_list_dep = []
    while 1:
        x = [x for x in G.nodes_iter() if G.in_degree(x)==0]
        #print '-->'," ".join(x)
        if x == []: break
        file_list_dep.extend(x)
        G.remove_nodes_from(x)
 
    # write makefile
    print "Writing makefile..."
    mf = open(os.path.join(srcOutDir,"makefile"),"w")
 
    # write options
    mf.write("BINDIR = %s\n"%bindir)
    mf.write("PROGRAM = %s\n"%program)
    mf.write("F90FLAGS = %s\n"%f90flags)
    mf.write("F90 = %s\n"%f90)
    if gperf:
        mf.write("USRLIBS = -lprofiler\n")
    mf.write("OBJECTS = \\\n")
   
    # remove header files
    tmp = copy.deepcopy(file_list_dep) 
    for f in tmp:
        if os.path.splitext(f)[1] not in srcExt:
            print "Removing %s.."%f
            file_list_dep.remove(f) 
    tmp = copy.deepcopy(file_list)  
    for f in tmp:
        if os.path.splitext(f)[1] not in srcExt:
            print "Removing %s.."%f
            file_list.remove(f)
 
    # first, write the depending files
    for f in file_list_dep:
        mf.write("\t%s%s \\\n"%(os.path.splitext(f)[0],'.o'))
        if f in file_list:
            file_list.remove(f)
    
    # second, write the independent files 
    for f in file_list[:-1]:
        mf.write("\t%s%s \\\n"%(os.path.splitext(f)[0],'.o'))
    f = file_list[-1]
    mf.write("\t%s%s \n"%(os.path.splitext(f)[0],'.o')) 

    # write tasks
    mf.write("all: %s\n"%program)
    mf.write("%s: $(OBJECTS)\n"%program)
    mf.write("\t-$(F90) $(F90FLAGS) -o $(BINDIR)$(PROGRAM) $(OBJECTS) $(USRLIBS) $(SYSLIBS)\n")
    for e in srcExt:
        mf.write("%%.o: %%%s\n"%e)
        mf.write("\t$(F90) $(F90FLAGS) -c $<\n")

    # write cleaning
    mf.write("clean:\n")
    mf.write("\trm -f *.o *.mod") 

    # close the make file
    mf.close()


