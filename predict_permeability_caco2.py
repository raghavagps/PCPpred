import argparse
import pandas as pd
import joblib
import os
import sys
import numpy as np
import pandas as pd
import numpy as np
import glob
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdMolDescriptors
from padelpy import padeldescriptor
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Descriptors3D
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import KFold
# from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, AdaBoostRegressor, ExtraTreesRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import LinearRegression
import lightgbm as lgb
import xgboost as xgb
import shutil

from collections import Counter
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=UserWarning)  
os.environ["LOKY_MAX_CPU_COUNT"] = "4"

def parse_args():
    parser = argparse.ArgumentParser(description="Predict permeability from SMILES input.")
    parser.add_argument("--input", type=str, required=True, help="Path to input file containing SMILES.")
    # parser.add_argument("--model", type=str, required=True, help="Meta model name for prediction.")
    parser.add_argument("--output", type=str, help="Path to save output CSV file")

    return parser.parse_args()

def load_smiles(input_path):
    try:
        with open(input_path, "r") as f:
            smiles = [line.strip() for line in f if line.strip()]
        return smiles
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

def create_smiles_csv(smiles_list):
    ids = [f"{i+1}" for i in range(len(smiles_list))]
    df = pd.DataFrame({
        'ID': ids,
        'SMILES': smiles_list
    })
    return df

# Mordred Descriptors
from mordred import Calculator, descriptors
calc_mord = Calculator(descriptors, ignore_3D=True)

def calculate_2dmordred_descriptors(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return calc_mord(mol)
        except Exception as e:
            print(f"Error processing SMILES {smiles} for Mordred descriptors: {e}")
            return None

def main():
    args = parse_args()
    smiles = load_smiles(args.input)
    # print(smiles)

    if not smiles:
        sys.exit("No SMILES found in input file.")
    
    test_df = create_smiles_csv(smiles)

    try:
        cwd = os.getcwd()

        descriptor_data = []
        for smiles in test_df['SMILES']:
            descriptors = calculate_2dmordred_descriptors(smiles)
            if descriptors is not None:
                descriptor_data.append(descriptors)
            else:
                descriptor_data.append([np.nan] * len(calc_mord.descriptors))
        
        raw_names = [desc.__class__.__name__ for desc in calc_mord.descriptors]
        counts = Counter()
        deduped_names = []
        for name in raw_names:
            if counts[name] == 0:
                deduped_names.append(name)
            else:
                deduped_names.append(f"{name}.{counts[name]}")
            counts[name] += 1

        # Convert Mordred results to DataFrame
        descriptor_df = pd.DataFrame(descriptor_data, columns=deduped_names)
        test_mordred_2d = pd.concat([test_df[['ID', 'SMILES']], descriptor_df], axis=1)
        # test_mordred_2d.to_csv('temp/test_2d_mordred.csv', index=False)

        
        models_dir = 'models/Caco2' 
        scaler_path = 'models/Caco2/scaler_caco2.joblib'                          
        model_base_name = 'LGBMRegressor' 
        features = joblib.load('models/Caco2/features.joblib')                  
        n_folds = 5  
         
        X_test = test_mordred_2d[features]
        obj_cols = X_test.select_dtypes(include=['object']).columns
        if len(obj_cols) > 0:
            for col in obj_cols:
                X_test[col] = 0

        
        scaler = joblib.load(scaler_path)
        X_new_scaled = scaler.transform(X_test)
        X_new_scaled = pd.DataFrame(X_new_scaled, columns= X_test.columns,index= X_test.index)
        all_fold_preds = []

        for fold in range(1, n_folds + 1):
            fold_model_path = os.path.join(models_dir, f"{model_base_name}_fold{fold}_Caco2.joblib")
            fold_model = joblib.load(fold_model_path)
            preds = fold_model.predict(X_new_scaled)
            preds = np.clip(preds, -10, -3.4)  
            all_fold_preds.append(preds)


        all_fold_preds = np.array(all_fold_preds)
        final_predictions = np.mean(all_fold_preds, axis=0)
        # print(final_predictions)
            
        output_df = pd.DataFrame({
        'SMILES': test_df['SMILES'],
        'Permeability': final_predictions
        })
        output_file = args.output if args.output else os.path.join(cwd, 'results', 'output_caco2.csv')
        
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
        output_df.to_csv(output_file, index=False)
    
        print("SMILES,Permeability")  
        for idx, row in output_df.iterrows():
            print(f"{row['SMILES']},{row['Permeability']:.2f}")
        
        print(f"Results saved to {output_file}")

    except Exception as e:
        sys.exit(f"Error during prediction: {str(e)}")
    
    finally:
        
        temp_dir = 'temp'
        if os.path.exists(temp_dir):
            try:
                for item in os.listdir(temp_dir):
                    item_path = os.path.join(temp_dir, item)
                    if os.path.isfile(item_path):
                        os.remove(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                
            except Exception as e:
                print(f"Warning: Could not clean up contents of {temp_dir}: {str(e)}")   


if __name__ == "__main__":
    main()
