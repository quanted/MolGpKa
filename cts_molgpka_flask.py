from flask import Flask, request
import logging
import json
import time

from cts_molgpka import CTSMolgpka


molgpka = CTSMolgpka()

app = Flask(__name__)


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("molgpka")
logger.info("API wrapper for molgpka.")


@app.route('/molgpka/test')
def status_check():
    return "molgpka running.", 200

@app.route('/molgpka/data')
def get_data():
    """
    Returns pka data from molgpka library.
    """
    t0 = time.time()
    args = request.args
    
    
    # TODO: User input data validation!


    if "smiles" in args:
        smiles = args["smiles"]
    else:
        return {"error": "Missing required molgpka parater 'smiles'"}, 400
    
    try:
        smiles, num_sites, pka_list, molgpka_smiles, pka_dict, molgpka_index = molgpka.main(smiles)
        results = {
            "status": True,
            "smiles": smiles,
            "num_sites": num_sites,
            "pka_list": pka_list,
            "pka_dict": pka_dict,
            "molgpka_smiles": molgpka_smiles,
            "molgpka_index": molgpka_index
        }

        # results = {"test": "testing"}
        # results = {'status': True, 'smiles': 'CC(O)=O', 'num_sites': 1, 'pka_list': [8.337571]}

        logging.warning("Pka list: {}".format(results["pka_list"]))
        logging.warning("Pka list: {}".format(type(results["pka_list"])))

        logging.warning("RESULTS: {}".format(results))

    except Exception as e:
        logging.error("molgpka_flask exception: {}".format(e))
        return {"error": "molgpka internal error"}, 500

    return results, 200



if __name__ == "__main__":
    logging.info("Starting up molgpka flask app")
    app.run(debug=True, host='0.0.0.0', port=8080)
