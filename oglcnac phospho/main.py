import argparse
import backend
import logging
import pandas as pd

def setup_logger(log_path: str, name: str = "logger") -> logging.Logger:
    """Configures the root logger.
    Parameters
    ----------
    log_path : str
        The filepath to the log handler.
    name : str (default: "logger")
        The name of the root logger.
    Returns
    -------
    logging.Logger
        The root logger.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename=log_path, encoding="utf-8", mode="w")
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def main():
    parser = argparse.ArgumentParser(prog = "command line tool")
    parser.add_argument("-n", "--nocache", action= "store_false", help= "use the cache files or re-pull from the API")
    parser.add_argument(
         "-o", "--output", type=str, default="tsv", nargs= "?", choices= ["csv", "tsv"], help="Output file type."
     )
    parser.add_argument(
         "-a", "--api", type=str, default="prod", nargs= "?", choices= ["beta", "prod"], help="API call method."
     )


    options= parser.parse_args()
    
    prod_flag = True 

    if options.api == "beta":
        prod_flag = False

    logger = setup_logger("./logfile.log")

    backend.get_all_list_data(options.nocache, logger, prod_flag)
    df = pd.read_csv("./output/supersearch_results.tsv", delimiter="\t")
    backend.get_protein_data(options.nocache, logger)
    #backend.get_protein_data(options.nocache, logger)
    backend.sites_data(logger, options.output, df)



if __name__ == "__main__":
    main()
