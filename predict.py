"""
    Predict labels for new samples
    Version 2.2.1
    Author:
        - Behnam Yousefi (behnm.yousefi@zmnh.uni-hamburg.de; yousefi.bme@gmail.com)
    Scientific Contributor:
        - Polina Gundorova (p.gundorova@uke.de)
    UKE, Hamburg, Germany - Jul 2024
"""

# Load the packages
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
import argparse

if __name__ == "__main__":

    # Read inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--train_data_dir', required=True, help='Train data directory; the output of the clustering pipeline provided in R')
    parser.add_argument('--new_data_dir', required=True, help='New data directory; the output of the gm.py for the train data')
    parser.add_argument('--k', default="3", help='Number of neighbours k in KNN.')

    # Set inputs
    args = parser.parse_args()
    train_data_dir = args.train_data_dir
    new_data_dir = args.new_data_dir
    k = args.k
    k = int(k)
    
    # Read the train data
    print("Step 1/6: Read the train data")
    df_train = pd.read_csv(train_data_dir)
    X = df_train[['Max_theory', 'Max_x', 'Max_y']]
    y = df_train['clusters']
    X, y = np.array(X), np.array(y)
    
    # Read the new data
    print("Step 2/6: Read the new data")
    df_test = pd.read_csv(new_data_dir)
    X_test = df_test[['Max_theory', 'Max_x', 'Max_y']]
    X_test = np.array(X_test)
    
    # Normalization
    print("Step 3/6: Normalization")
    m, s = X.mean(axis=0), X.std(axis=0)
    X = (X - m) / s
    X_test = (X_test - m) / s
    
    # Set (train) the model
    print("Step 4/6: Train the model")
    neigh = KNeighborsClassifier(n_neighbors=k)
    neigh.fit(X, y)
    
    # Infere the output for new samples (test)
    print("Step 5/6: Infere the output for new samples")
    y_hat = neigh.predict(X_test)
    
    # Save the predicted labels
    print("Step 6/6: Save the predicted labels")
    df_test['predicted_lanbel'] = y_hat
    df_test.to_csv(new_data_dir.split('.csv')[0] + '_predict.csv')
    print("The output is saved at:")
    print(new_data_dir.split('.csv')[0] + '_predict.csv')
