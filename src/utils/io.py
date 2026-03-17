import json
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    """ json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.int64, np.int32, np.int8)):
            return int(obj)
        if isinstance(obj, (np.float64, np.float32)):
            return float(obj)
        return json.JSONEncoder.default(self, obj)

def save_json(data, path):
    with open(path, 'w') as f:
        json.dump(data, f, cls=NumpyEncoder)