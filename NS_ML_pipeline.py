#!/usr/bin/env python
#Author: Behzad Moumbeini


#import required packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.cluster import KMeans
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
from keras.models import Sequential
from keras.layers import Dense
from keras.callbacks import EarlyStopping



#fraction is the ration for test set
def make_train_test(infile, fraction):

    data = pd.read_csv(infile, sep="\t")
    keep = ['CodeClass','Name', 'Accession']

    columns = data.columns.drop(keep)
    cut = max(1, int((1-fraction)*len(columns)))

    train = columns[:cut] # unique headers for train
    test = columns[cut:]  # unique headers for test

    full_train = data[ keep + train.to_list() ]
    full_test = data[ keep + test.to_list() ]
    full_train.to_csv('full_train.csv', header=True ) 
    full_test.to_csv('full_test.csv', header=True ) 



def process_metada(metadata_file):
    y = pd.read_csv(metadata_file, sep="\t")
    y['group'].value_counts()
    labels = ['A', 'B']
    return [labels, y]



def process(full_train, full_test):
    df_train = pd.read_csv(full_train, sep=',')
    df_test = pd.read_csv(full_test, sep=",")

    train_to_keep = [col for col in df_train.columns if "call" not in col]
    test_to_keep = [col for col in df_test.columns if "call" not in col]
    X_train_tr = df_train[train_to_keep]
    X_test_tr = df_test[test_to_keep]
    return [X_train_tr, X_test_tr]


def transpose(X_train_tr, X_test_tr):
    X_train = X_train_tr.T
    X_test = X_test_tr.T

    X_train.columns = X_train.iloc[1]
    X_train = X_train.drop(['Unnamed: 0', 'CodeClass', 'Name', 'Accession']).apply(pd.to_numeric)

    # Clean up the column names for Testing data
    X_test.columns = X_test.iloc[1]
    X_test = X_test.drop(['Unnamed: 0', 'CodeClass', 'Name', 'Accession']).apply(pd.to_numeric)
    return [X_train, X_test]



def pca3(X_train):
    pca3 = PCA(n_components=3).fit(X_train)
    X_train_reduced = pca3.transform(X_train)
    return X_train_reduced

def add_column(X_train, X_test):
    p_train = X_train.index.values
    X_train.insert( 0, column="CodeClass2",value = p_train)
    p_test = X_test.index.values
    X_test.insert( 0, column="CodeClass2",value = p_test)
    return [X_train2, X_test2]


def make_y_train_test(X_train, X_test, y):
    y_train = y[y['sample'].isin(X_train['CodeClass2'])]
    y_test = y[y['sample'].isin(X_test['CodeClass2'])]
    return [y_train, y_test]


def scatter(X_train_reduced):
    color = y_train.iloc[:, 1].replace({'B': 2, 'A': 1})
    fig = plt.figure(1, figsize = (10, 6))
    plt.figure()
    plt.scatter(X_train_reduced[:, 0],  X_train_reduced[:, 1] , color, cmap = plt.cm.Paired , linewidths=10)
    plt.annotate('Note the Brown Cluster', xy = (30000,-2000))
    plt.title("2D Transformation of the Above Graph ")
    plt.savefig('books_read.pdf')



def PCA_3D(X_train_reduced):
    plt.clf()
    fig = plt.figure(1, figsize=(10,6 ))
    ax = Axes3D(fig, elev=-150, azim=110,)
    ax.scatter(X_train_reduced[:, 0], X_train_reduced[:, 1], X_train_reduced[:, 2], cmap = plt.cm.Paired, linewidths=10)
    ax.set_title("First three PCA directions")
    ax.set_xlabel("1st eigenvector")
    ax.w_xaxis.set_ticklabels([])
    ax.set_ylabel("2nd eigenvector")
    ax.w_yaxis.set_ticklabels([])
    ax.set_zlabel("3rd eigenvector")
    ax.w_zaxis.set_ticklabels([])
    plt.savefig('read.pdf')



def kmeans(X_train_tr, X_test_tr):
    X_train = transpose(X_train_tr, X_test_tr)[0]
    X_train_fl = X_train.astype(float, 64)
    X_test = transpose(X_train_tr, X_test_tr)[1]
    X_test_fl = X_test.astype(float, 64)

    # Apply the same scaling to both datasets
    scaler = StandardScaler()
    X_train_scl = scaler.fit_transform(X_train_fl)
    X_test_scl = scaler.transform(X_test_fl)
    kmeans = KMeans(n_clusters=2, random_state=0).fit(X_train_scl)
    km_pred = kmeans.predict(X_test_scl)
    accuracy = 'K-means accuracy:', round(accuracy_score(y_test.iloc[:, 1].replace({'B': 0, 'A': 1}), km_pred), 3)
    confusion_matrix_km = confusion_matrix(y_test.iloc[:, 1].replace({'B': 0, 'A': 1}), km_pred)

    return [km_pred, accuracy, confusion_matrix_km, X_train_scl, X_test_scl]




def kmcm(confusion_matrix_km):

    ax = plt.subplot()
    sns.heatmap(confusion_matrix_km, annot=True, ax = ax, fmt='g', cmap='Greens') 

    # labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels') 
    ax.set_title('K-means Confusion Matrix') 
    ax.xaxis.set_ticklabels(labels) 
    ax.yaxis.set_ticklabels(labels, rotation=360);
    ax.figure.savefig("confusion_matrix.pdf")

# Create a Gaussian classifier
def Gaussian_classifier(X_train,y_train, X_test, y_test):
    
    nb_model = GaussianNB()

    nb_model.fit(X_train.iloc[:,2:], y_train.iloc[:,1])

    nb_pred = nb_model.predict(X_test.iloc[:,2:])

    print('Naive Bayes accuracy:', round(accuracy_score(y_test.iloc[:, 1].replace({'B': 0, 'A': 1}), nb_pred), 3))

    cm_nb =  confusion_matrix(y_test.iloc[:,1], nb_pred)

    ax = plt.subplot()
    sns.heatmap(cm_nb, annot=True, ax = ax, fmt='g', cmap='Greens') 

    # labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels') 
    ax.set_title('Naive Bayes Confusion Matrix') 
    ax.xaxis.set_ticklabels(labels) 
    ax.yaxis.set_ticklabels(labels, rotation=360);
    ax.figure.savefig("confusion_matrix_Gaussian.pdf")



def logistic_regression(X_train,y_train, X_test, y_test):

    log_grid = {'C': [1e-03, 1e-2, 1e-1, 1, 10], 
                     'penalty': ['l1', 'l2']}

    log_estimator = LogisticRegression(solver='liblinear')

    log_model = GridSearchCV(estimator=log_estimator, 
                      param_grid=log_grid, 
                      cv=3,
                      scoring='accuracy')

    log_model.fit(X_train.iloc[:,2:], y_train.iloc[:,1])

    print("Best Parameters:\n", log_model.best_params_)

    # Select best log model
    best_log = log_model.best_estimator_

    # Make predictions using the optimised parameters
    log_pred = best_log.predict(X_test.iloc[:,2:])

    print('Logistic Regression accuracy:', round(accuracy_score(y_test.iloc[:,1], log_pred), 3))

    cm_log =  confusion_matrix(y_test.iloc[:,1], log_pred)

    ax = plt.subplot()
    sns.heatmap(cm_log, annot=True, ax = ax, fmt='g', cmap='Greens') 

    # labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels') 
    ax.set_title('Logistic Regression Confusion Matrix') 
    ax.xaxis.set_ticklabels(labels) 
    ax.yaxis.set_ticklabels(labels, rotation=360);
    ax.figure.savefig("confusion_matrix_logistic_regression.pdf")


def SVM(X_train_pca,y_train, X_test_pca, y_test):


    # Parameter grid
    svm_param_grid = {'C': [0.1, 1, 10, 100], 'gamma': [1, 0.1, 0.01, 0.001, 0.00001, 10], "kernel": ["linear", "rbf", "poly"], "decision_function_shape" : ["ovo", "ovr"]} 

    # Create SVM grid search classifier
    svm_grid = GridSearchCV(SVC(), svm_param_grid, cv=3)

    # Train the classifier
    svm_grid.fit(X_train_pca, y_train.iloc[:,1])

    print("Best Parameters:\n", svm_grid.best_params_)

    # Select best svc
    best_svc = svm_grid.best_estimator_

    # Make predictions using the optimised parameters
    svm_pred = best_svc.predict(X_test_pca)

    print('SVM accuracy:', round(accuracy_score(y_test.loc[:,1], svm_pred), 3))

    cm_svm =  confusion_matrix(y_test.iloc[:,1], svm_pred)

    ax = plt.subplot()
    sns.heatmap(cm_svm, annot=True, ax = ax, fmt='g', cmap='Greens') 

    # Labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels') 
    ax.set_title('SVM Confusion Matrix') 
    ax.xaxis.set_ticklabels(labels) 
    ax.yaxis.set_ticklabels(labels, rotation=360);
    ax.figure.savefig("confusion_matrix_SVM.pdf")


def random_forest(X_train, y_train, X_test, y_test):

    # Hyperparameters search grid 
    rf_param_grid = {'bootstrap': [False, True],
             'n_estimators': [60, 70, 80, 90, 100],
             'max_features': [0.6, 0.65, 0.7, 0.75, 0.8],
             'min_samples_leaf': [8, 10, 12, 14],
             'min_samples_split': [3, 5, 7]
            }

    # Instantiate random forest classifier
    rf_estimator = RandomForestClassifier(random_state=0)

    # Create the GridSearchCV object
    rf_model = GridSearchCV(estimator=rf_estimator, param_grid=rf_param_grid, cv=3, scoring='accuracy')

    # Fine-tune the hyperparameters
    rf_model.fit(X_train.iloc[:,2:], y_train.iloc[:,1])

    print("Best Parameters:\n", rf_model.best_params_)

    # Get the best model
    rf_model_best = rf_model.best_estimator_

    # Make predictions using the optimised parameters
    rf_pred = rf_model_best.predict(X_test.iloc[:,2:])

    print('Random Forest accuracy:', round(accuracy_score(y_test.iloc[:,1], rf_pred), 3))

    cm_rf = confusion_matrix(y_test.iloc[:,1], rf_pred)

    ax = plt.subplot()
    sns.heatmap(cm_rf, annot=True, ax = ax, fmt='g', cmap='Greens') 

    # labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels') 
    ax.set_title('Random Forest Confusion Matrix') 
    ax.xaxis.set_ticklabels(labels) 
    ax.yaxis.set_ticklabels(labels, rotation=360);
    ax.figure.savefig("confusion_matrix_random_forest.pdf")



# XGB — PCA with Grid Search
def XGB_PCA(X_train_pca, y_train, X_test_pca, y_test):


    xgb_grid_params = {'max_depth': [3, 4, 5, 6, 7, 8, 10, 12],
                   'min_child_weight': [1, 2, 4, 6, 8, 10, 12, 15],
                   'n_estimators': [40, 50, 60, 70, 80, 90, 100, 110, 120, 130],
                   'learning_rate': [0.001, 0.01, 0.05, 0.1, 0.2, 0.3]}

    fixed_params = {'random_state': 0,
                    'n_jobs': -1}

    xgb_model = GridSearchCV(xgb.XGBClassifier(**fixed_params), 
                           param_grid = xgb_grid_params, 
                           scoring = 'accuracy',
                           cv = 3)

    xgb_model.fit(X_train_pca, y_train.iloc[:,1])

    print("Best Parameters:\n", xgb_model.best_params_)

    # Get the best model
    xgb_model_best = xgb_model.best_estimator_

    # Make predictions using the optimised parameters
    xgb_pred = xgb_model_best.predict(X_test_pca)

    print('XGB (PCA with Grid Search) accuracy:', round(accuracy_score(y_test.iloc[:,1], xgb_pred), 3))

    cm_xgb = confusion_matrix(y_test.iloc[:,1], xgb_pred)

    ax = plt.subplot()
    sns.heatmap(cm_xgb, annot=True, ax = ax, fmt='g', cmap='Greens') 

    # labels, title and ticks
    ax.set_xlabel('Predicted labels')
    ax.set_ylabel('True labels') 
    ax.set_title('XGB (PCA with Grid Search) Confusion Matrix') 
    ax.xaxis.set_ticklabels(labels) 
    ax.yaxis.set_ticklabels(labels, rotation=360);
    ax.figure.savefig("confusion_matrix_XGB.pdf")



if __name__ == "__main__":

        separation = make_train_test("Normalized_data.txt", 0.3)
        labels = process_metada("metadata.txt")[0]
        y = process_metada("metadata.txt")[1]
        X_train_tr = process("full_train.csv", "full_test.csv")[0]
        X_test_tr = process("full_train.csv", "full_test.csv")[1]
        X_train = transpose(X_train_tr, X_test_tr)[0]
        X_test = transpose(X_train_tr, X_test_tr)[1]
        pca = PCA()
        pca.fit_transform(X_train)


        total = sum(pca.explained_variance_)
        k = 0
        current_variance = 0
        while current_variance/total < 0.90:
         current_variance += pca.explained_variance_[k]
         k = k + 1

        print(k, " features explain around 90% of the variance. From 7129 features to ", k, ", not too bad.", sep='')

        pca = PCA(n_components=k)
        X_train.pca = pca.fit(X_train)
        X_train_pca = pca.transform(X_train)
        X_test_pca = pca.transform(X_test)

        var_exp = pca.explained_variance_ratio_.cumsum()
        var_exp = var_exp*100
        plt.bar(range(k), var_exp);

        X_train_reduced = pca3(X_train)
        X_train2 = add_column(X_train, X_test)[0]
        X_test2 = add_column(X_train, X_test)[1]

        y_train = make_y_train_test(X_train, X_test, y)[0]
        y_test = make_y_train_test(X_train, X_test, y)[1]

        print(scatter(X_train_reduced))
        print(PCA_3D(X_train_reduced))
        print(type(y_test.iloc[:,1]))
        y2 = y_test.iloc[:, 1].replace({'B': 0, 'A': 1})

        km_pred = kmeans(X_train_tr, X_test_tr)[0]
        accuracy = kmeans(X_train_tr, X_test_tr)[1]
        confusion_matrix_km = kmeans(X_train_tr, X_test_tr)[2]
        X_train_scl = kmeans(X_train_tr, X_test_tr)[3]
        X_test_scl = kmeans(X_train_tr, X_test_tr)[4]
        print(kmcm(confusion_matrix_km))
        print(Gaussian_classifier(X_train, y_train, X_test, y_test))
        print(logistic_regression(X_train,y_train, X_test, y_test))
        print(SVM(X_train_pca,y_train, X_test_pca, y_test))
        print(random_forest(X_train, y_train, X_test, y_test))
        print(XGB_PCA(X_train_pca, y_train, X_test_pca, y_test))

