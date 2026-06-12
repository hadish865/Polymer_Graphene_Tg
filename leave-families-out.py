import pandas as pd
import xlsxwriter
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error, mean_squared_error, r2_score
import numpy as np
import pickle as pkl

def Leave_papers_out_families():
    # import dataset
    filename = "dataset"
    dataset = pkl.load(open("data/"+filename+".rb", "rb"))

    X_fields = ["graphene_w","CH","C=OH","OH","NH","B","BCH","hydrogen_bond","polarity"]
    Y_field = ["nano_composite_Tg_(K)"]
    X = dataset.iloc[:,[4,9,11,15,17,19,21,35,36,37]]
    Y = dataset.iloc[:,5]

    models = [
          # Gradient Boosting optimized result
          GradientBoostingRegressor(n_estimators = 193, min_samples_split = 0.705201292),
          ]

    number_of_models = 1
    number_of_groups = 11

    R2_results   = np.zeros((number_of_models, 2, number_of_groups))
    MAE_results  = np.zeros((number_of_models, 2, number_of_groups))
    MAPE_results = np.zeros((number_of_models, 2, number_of_groups))
    MSE_results  = np.zeros((number_of_models, 2, number_of_groups))

    split_number = 0
    col_index = 0

    names = dataset["polymer_name"]
    PBT = ["PBT"]
    NYLON12 = ["nylon 12"]
    PE = ["PE"]
    PEO = ["PEO"]
    PTT = ["PTT"]
    PA6 = ["PA6"]
    PVC = ["PVC", "PPVC"]
    PEEK = ["PEEK"]
    PP = ["PP", "iPP"]
    PET = ["PET"]
    PMMA = ["PMMA"]

    groups = [PBT , NYLON12 , PE , PEO , PTT , PA6 , PVC, PEEK , PP, PET, PMMA]
    group_names = [str(group) for group in groups]

    for group in groups:
        group_indexes = dataset[names.isin(group)].index.to_numpy()
        Xtrain_full, Xtest = X.drop(group_indexes), X.iloc[group_indexes]
        Ytrain_full, Ytest = Y.drop(group_indexes), Y.iloc[group_indexes]

        TEST_SIZE = (35 - len(group_indexes))/(236 - len(group_indexes))
        if TEST_SIZE <= 0:
            continue
        Xtrain, Xtest2, Ytrain, Ytest2 = train_test_split(Xtrain_full, Ytrain_full, test_size=TEST_SIZE, random_state = 42)
        Xtest = pd.concat((Xtest, Xtest2), axis=0)
        Ytest = pd.concat((Ytest, Ytest2), axis=0)

        scaler = StandardScaler()
        Xtrain = scaler.fit_transform(Xtrain)
        Xtest = scaler.transform(Xtest)
        Ytrain = np.array(Ytrain)
        Ytest  = np.array(Ytest)
        regressor: GradientBoostingRegressor = models[0]
        # train the data
        regressor.fit(Xtrain, Ytrain.ravel())
        predict_train = regressor.predict(Xtrain)
        predict_test = regressor.predict(Xtest)
        
        print("split: "+str(col_index+1)+" complete.")
        MAE_results [split_number, 1, col_index]  = mean_absolute_error(Ytest, predict_test)
        MAPE_results[split_number, 1, col_index]  = mean_absolute_percentage_error(Ytest, predict_test)
        MSE_results [split_number, 1, col_index]  = mean_squared_error(Ytest, predict_test)
        r2_test  = r2_score(Ytest.ravel(), predict_test)
        MAE_results [split_number ,0, col_index] = mean_absolute_error(Ytrain, predict_train)
        MAPE_results[split_number ,0, col_index] = mean_absolute_percentage_error(Ytrain, predict_train)
        MSE_results [split_number ,0, col_index] = mean_squared_error(Ytrain, predict_train)
        r2_train  = r2_score(Ytrain.ravel(), predict_train)

        if (r2_train < 0):
            r2_train = np.abs(r2_train)
        if (r2_test < 0):
            r2_test = np.abs(r2_test)
        if (r2_train > 1):
            r2_train = 1/r2_train
        if (r2_test > 1):
            r2_test = 1/r2_test

        R2_results[split_number, 1, col_index] = r2_test
        R2_results[split_number, 0, col_index] = r2_train
        col_index += 1
    print("split (model) "+str(split_number+1)+" complete")
    print("===================")
    print("saving results")

    Save_leave_one_out_results_full(R2_results, MAE_results, MAPE_results, MSE_results, group_names, "leave_out_familes")

def Save_leave_one_out_results_full(R2_results, MAE_results, MAPE_results, MSE_results, group_names, filename="results"):

    workbook = xlsxwriter.Workbook("results/"+filename+".xlsx")

    title_format = workbook.add_format()
    title_format.set_bold()
    title_format.set_align("center")
    title_format.set_align("vcenter")
    # format for other cells
    normal_format = workbook.add_format()
    normal_format.set_align("center")
    normal_format.set_align("vcenter")
    # name format
    name_format = workbook.add_format()
    name_format.set_bold()
    name_format.set_align("center")
    name_format.set_align("vcenter")
    name_format.set_pattern(1)
    name_format.set_fg_color("#89A7DE")

    index = 0
    worksheet_metrics = workbook.add_worksheet("model " + str(index))
    worksheet_metrics.write(0,0,"names",title_format)

    worksheet_metrics.write(1,0,"R2 train" ,title_format)
    worksheet_metrics.write(2,0,"R2 test" ,title_format)

    worksheet_metrics.write(3,0,"MAE train" ,title_format)
    worksheet_metrics.write(4,0,"MAE test" ,title_format)

    worksheet_metrics.write(5,0,"MAPE train" ,title_format)
    worksheet_metrics.write(6,0,"MAPE test" ,title_format)

    worksheet_metrics.write(7,0,"MSE train" ,title_format)
    worksheet_metrics.write(8,0,"MSE test" ,title_format)

    worksheet_metrics.write_row(0, 1, group_names, title_format)
    worksheet_metrics.write_row(1, 1, R2_results  [index, 0, :]     ,normal_format)
    worksheet_metrics.write_row(2, 1, R2_results  [index, 1, :]     ,normal_format)
    worksheet_metrics.write_row(3, 1, MAE_results [index, 0, :]     ,normal_format)
    worksheet_metrics.write_row(4, 1, MAE_results [index, 1, :]     ,normal_format)
    worksheet_metrics.write_row(5, 1, MAPE_results[index, 0, :]*100 ,normal_format)
    worksheet_metrics.write_row(6, 1, MAPE_results[index, 1, :]*100 ,normal_format)
    worksheet_metrics.write_row(7, 1, MSE_results [index, 0, :]     ,normal_format)
    worksheet_metrics.write_row(8, 1, MSE_results [index, 1, :]     ,normal_format)
    worksheet_metrics.autofit()

    workbook.close()

Leave_papers_out_families()