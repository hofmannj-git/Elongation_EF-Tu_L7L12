import numpy as np
import sys
import h5py as h5
import tqdm
import os

# organize boolean arrays by time
files = np.array(os.listdir('search'))
times = np.array([int(x.split('.')[0][7:]) for x in files])
sort_idx = np.argsort(times)
times = times[sort_idx]
files = files[sort_idx]

# Create functions
def find_unique_res(array):
    """
    Find unique values and list of indices corresponding to repeats of that value
    
    INPUT
    -----
    array = np sequence of values
    
    
    OUTPUT
    ------
    vals = np.array of unique values of array
    res = len(vals) list of indices corresponding to unique vallues
    """

    # creates an array of indices, sorted by unique element
    idx_sort = np.argsort(array)

    # sorts records array so all unique elements are together 
    sorted_records_array = array[idx_sort]

    # returns the unique values, the index of the first occurrence of a value, and the count for each element
    uni, idx_start, count = np.unique(sorted_records_array, return_counts=True, return_index=True)

    # splits the indices into separate arrays
    uni_res = np.split(idx_sort, idx_start[1:])
    
    return uni, uni_res,idx_start

def consecutive(data, stepsize=1):
    """
    Find values that are in consectuive order in an array
    
    INPUT
    -----
    array = np sequence of values
    stepsize = optional distance between steps
    
    OUTPUT
    ------
    consec = list of values from data that are consecutive 
    where = np array telling where to split array to find new consecutive sequences
    """
    return np.split(data,np.where(np.diff(data) != stepsize)[0]+1),(np.where(np.diff(data) != stepsize)[0]+1).tolist()

def rebind(r1,r2):
    """
    Helper function to find if two arrays have any intersects
    """
    if len(np.intersect1d(r1,r2)) > 0:
        return True
    else:
        return False

def list_splitter(l,split_idx):
    if len(split_idx) == 0:
        return [l]
    
    else:
        lol = []
        split_idx_buffered = np.concatenate([np.array([0]),split_idx,np.array([len(l)])])
        for i in range(len(split_idx_buffered)-1):
            lol.append(l[split_idx_buffered[i]:split_idx_buffered[i+1]])
        return lol
def get_bind_search(t_id,r_id):
    """
    Function for getting the binding and search times of pairwise interactions. Immediate Rebinding is treated as discrete 
    event. Search times ignore immedaite repeat binding.
    
    INPUT
    -----
    
    t_id = np sequence of times where interactions occur
    r_id = np sequence of indices of receptors
    
    OUTPUT
    ------
    search = int list of search times, ignoring immediate repeat binding
    l_bind = int list of binding times, treating immediate repeat binding as a discrete binding event
    partners = list of list of binding partners
    t_split = list of list of times grouped by binding event.
    """
    
    # find all unique times in array
    uni, uni_res,idx_start = find_unique_res(t_id)

    # separate binding events by where jumps occur
    _,split_idx = consecutive(uni)

    # find unbinding events where a time jump does not occur
    rib_uni = [r_id[x].tolist() for x in uni_res]
    unbind_z = np.where(np.array([len(np.intersect1d(x,rib_uni[i+1])) for i, x in enumerate(rib_uni[:-1])]) == 0)[0]+1
    for x in unbind_z.tolist():
        split_idx.extend([x])

    split_idx = np.sort(np.unique(np.array(split_idx)))
    t_split =  np.split(uni,split_idx)
    
    # lengths of binding times
    l_bind = [len(b) for b in t_split]

    # lengths of unbound times
    unbind = [t_split[i+1][0] - x[-1]-1 for i,x in enumerate(t_split[:-1])] 

    # get receptor id for each binding event
    partners = [r_id[x].tolist() for x in uni_res]
    partners_split = list_splitter(partners,split_idx)
    
    # filter for search times between non-identical receptors
    search = []
    for i,x in enumerate(unbind):
        if rebind(partners_split[i][-1],partners_split[i+1][0]):
            continue
        else:
            search.append(unbind[i])
            
    return search, l_bind, partners,t_split


# define number of snap interval (tlen) to read per h5 file, and number of iterations (iters) to get through the h5 file 
n_tc = 7476
tlen = 100
iters = 1
bound_all = [[] for x in range(n_tc)]
search_all = [[] for x in range(n_tc)]
files_to_read = files[-80:]

# Define number of TCs

# define number of snap interval (tlen) to read per h5 file, and number of iterations (iters) to get through the h5 file 
tlen = 10
iters = 1
bound_all = [[] for x in range(n_tc)]
search_all = [[] for x in range(n_tc)]
tsplits = [[] for idx in range(n_tc)]
css = [True for idx in range(n_tc)]
pts = [[] for idx in range(n_tc)]
for t,f in enumerate(files_to_read):
    print(f)
    bound_idx = np.load('search/' +f)[:-1]
    
    for i in range(n_tc):
        search1 = []
        
            
        # find where binding events occur
        t_id,r_id = np.where(bound_idx[:tlen,:,i])
        
        # decide if first time is search event
        if t == 0 and len(t_id) > 0 and t_id[0] == 0 and css[i]:
            css[i] = False
        
        if t > 0:
            
            # rectify time index to account for reindexing each iteration
            t_id += tlen*t
            
            # add back last binding event from last iteration if it exists
            if len(tsplits[i]) > 0:

                tstart = []
                count = 0
                for j in tsplits[i]:
                    tstart.append(np.arange(count,count+len(j)))
                    count += len(j)
                r_add = [pts[i][y] for y in tstart[-1]]
                t_add = [[x]*len(r_add[fb]) for fb,x in enumerate(tsplits[i][-1])]
                t_add = [item for sublist in t_add for item in sublist]
                r_add = [item for sublist in r_add for item in sublist]
                t_id = np.concatenate([t_add,t_id])
                r_id = np.concatenate([r_add,r_id])

        # find binding partners
        search, l_bind, pts[i],tsplits[i] = get_bind_search(t_id,r_id)

        # find first search length as the time of first binding event
        if len(t_id) > 0 and t_id[0] > 0  and css[i]:
            search = [t_id[0]] + search
            css[i] = False
        
        # output bind and search values
        if t == len(files_to_read)-1:
            bound_all[i].append(l_bind)          
        else:
            bound_all[i].append(l_bind[:-1])
        search_all[i].append(search)

np.save('bound_t.npy',bound_all)
np.save('search_t.npy',search_all)



