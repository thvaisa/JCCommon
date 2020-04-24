import json



def read_json(fname,txt=None):
    if(txt is not None): print(txt.format(fname))
    with open(fname,"r") as json_file:
        return json.load(json_file)
    return None

def save_json(fname,task,txt=None):
    if(txt is not None): print(txt.format(fname))
    with open(fname, 'w') as outfile:
        json.dump(task, outfile)


