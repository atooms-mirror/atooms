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
