class NpEncoder(json.JSONEncoder):
    # encoder used to write json
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def export_uorfs(uorfs):
    """Write the output uorf JSON

    :param uorfs: JSON structure containing uORF data
    :type uorfs: JSON
    """
    # Ensure the directory exists
    os.makedirs('./data', exist_ok=True)
    
    # Write to JSON file
    with open('./data/uorfs.json', 'w') as f:
        #json.dump(uorfs, f)
        json.dump(uorfs, f, cls=NpEncoder)
    
