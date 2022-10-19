import requests
import json

def get_url_dockign_seamdock(path_mol_pdb, path_prot_pdb, website_param = True):
    """
    Create a new seamless instance and launch docking for a specified molecule and protein
    """
    url = 'https://seamless.rpbs.univ-paris-diderot.fr/cloudless/launch'
    rq = requests.post(url, {"service": "seamdock"}, allow_redirects=False)

    # ADD molecule
    url = 'https' + rq.headers["Location"][4:] + "ctx/data/ligand_input"
    url_type = 'https' + rq.headers["Location"][4:] + "ctx/data/ligand_ext.txt"

    mol_pdb = open(path_mol_pdb).read()
    headers = {
        "Content-Type": "application/json; charset=utf-8",
    }

    payload_type = json.dumps({"buffer": "pdb", "marker": 1})
    payload = json.dumps({"buffer": mol_pdb, "marker": 2})

    requests.put(url_type,
                 data=payload_type, headers=headers)
    requests.put(url,
                 data=payload, headers=headers)

    # ADD PDB
    url = 'https' + rq.headers["Location"][4:] + "ctx/data/receptor.pdb"
    url_type = 'https' + rq.headers["Location"][4:] + "ctx/data/rec_ext.txt"

    headers = {
        "Content-Type": "application/json; charset=utf-8",
    }

    pdb = open(path_prot_pdb).read()
    payload = json.dumps({"buffer": pdb, "marker": 3})
    payload_type = json.dumps({"buffer": "pdb", "marker": 4})

    requests.put(url_type,
                 data=payload_type, headers=headers)
    requests.put(url,
                 data=payload, headers=headers)

    # ADD BOC PARAMETERS
    if website_param :
        url = 'https' + rq.headers["Location"][4:] + "ctx/data/rec_view_opt.json"
        rec_view_opt_json = "{'box_alpha': 50, 'boc_color': '#000000', 'cartoon': false, 'licorice': true, 'surf_alpha': 0.5, 'surface': true, 'wireframe': false, 'x': 0, 'x_size': 43, 'y': 0, 'y_size': 24, 'z': 0, 'z_size': 38}"
        requests.put(url, data=json.dumps({"buffer": rec_view_opt_json, "marker": 30}), headers=headers)

    # GET MAIN URL
    url_seamdock = 'https' + rq.headers["Location"][4:] + "ctx/index.html"
    return(url_seamdock)