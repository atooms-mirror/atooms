def get_step_from_path(f):
    s = re.search(r'%s-(\d*)' % basename, f)
    if s:
        return int(s.group(1))

def sort_files_steps(files, steps):
    file_steps = zip(files, steps)
    file_steps.sort(key = lambda a : a[1])
    new_files = [a[0] for a in file_steps] 
    new_steps = [a[1] for a in file_steps]
    return new_files, new_steps

def get_period(data):
    if len(data) < 2:
        return 1
    delta_old = 0
    delta_one = data[1] - data[0]
    iold = data[0]
    period = 0
    for ii in range(1, len(data)):
        i = data[ii]
        delta = i-iold
        if delta < delta_old and delta == delta_one and abs(delta-delta_old)>2:
            return period
        else:
            period += 1
            iold = i
            delta_old = delta
    return 1
