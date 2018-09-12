# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:34:58 2018

@author: 77170
"""

import csv
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import levene
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.dummy import DummyClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import StratifiedKFold
import pandas as pd
from pandas import DataFrame
from sklearn.metrics import average_precision_score,precision_score,recall_score,f1_score
from sklearn.metrics import roc_auc_score,accuracy_score,matthews_corrcoef,roc_curve
import math
import matplotlib.pyplot as plt
import seaborn as sns
#将原数据中的cisplatin的数据提取出来写入pre_pati_data中
'''
pati_Data=open("D:/py_data/pati_data.csv",'rU')
csvReader=csv.reader(pati_Data)
sheet=[]
for row in csvReader:
    sheet.append(row)
pati_Data.close()
sheet_1=[]
for i in range(len(sheet)):
    if sheet[i][2]=="Cisplatin":
        sheet_1.append(sheet[i])
pre_pati_Data=open("D:/py_data/pre_pati_data.csv","wb")
csvWriter=csv.writer(pre_pati_Data)
for row in sheet_1:
    csvWriter.writerow(row)
    print row
pre_pati_Data.close()
'''
#将pre_pati_data中数据改成和miRNA数据一样的格式并且标注sensitive和insensitive
'''
pre_pati_Data=open('D:/py_data/pre_pati_data.csv','rU')
csvReader=csv.reader(pre_pati_Data)
sheet=[]
for row in csvReader:
    sheet.append(row)
pre_pati_Data.close()
Data=open('D:/py_data/data.csv','wb')
csvWriter=csv.writer(Data)
sheet_1=[]
sheet_1.append(['SampleID','Response'])
print sheet[0]
for i in range(len(sheet)):
    sheet[i][0]=sheet[i][0][sheet[i][0].index('(')+1:sheet[i][0].index(')')]+'.'+sheet[i][1].replace('-','.')+'.01'
    if sheet[i][4]=='Stable Disease' or sheet[i][4]=='Clinical Progressive Disease':
        sheet[i][1]='insensitive'
    else:
        sheet[i][1]='sensitive'
    sheet_1.append(sheet[i][:2])
for row in sheet_1:
    csvWriter.writerow(row)
Data.close()  
''' 
'''
#将新的TCGA的数据进行提取，处理和分开
new_data = open('D:/clinical/TCGA-BLCA/BLCA.csv','rU')
csvReader = csv.reader(new_data)
sheet = []
for row in csvReader:
    sheet_1 = []
    if row[3] == 'bcr_patient_barcode':
        sheet.append(['SampleID','drug','response'])
    else:
        sheet_1.append(row[3]+'-01')
        sheet_1.append(row[25].split('##'))
        sheet_1.append(row[26].split('##'))
        for i in range(len(sheet_1[1])):
            sheet_2 = []
            sheet_2.append(sheet_1[0])
            sheet_2.append(sheet_1[1][i])
            sheet_2.append(sheet_1[2][i])
            sheet.append(sheet_2)
new_data.close()
new_Data = open('D:/clinical/TCGA-BLCA/BLCA1.csv','wb')
csvWriter = csv.writer(new_Data)
for row in sheet:
    csvWriter.writerow(row)
new_Data.close()
'''
'''
#将处理后的数据提取我们需要的cisplatin数据
data = open('D:/clinical/TCGA-BLCA/BLCA1.csv','rU')
csvReader = csv.reader(data)
sheet = []
sheet.append(['SampleID','Response'])
for row in csvReader:
    sheet_1 = []
    if row[1] == 'cisplatin' or row[1] == 'Cisplatin':
        if row[2] == 'Complete Response' or row[2] == 'Partial Response':
            sheet_1.append(row[0])
            sheet_1.append('Sensitive')
        else:
            sheet_1.append(row[0])
            sheet_1.append('Insensitive')
    if sheet_1:
        sheet.append(sheet_1)
    else:
        pass
final_data = open('D:/clinical/TCGA-BLCA/BLCA2.csv','wb')
csvWriter = csv.writer(final_data)
for row in sheet:
    csvWriter.writerow(row)
final_data.close()
'''
'''
#将miRNA数据处理成需要的数据
f = open('D:/clinical/TCGA-BLCA/BLCA-miRNA/Merge_matrix.csv','rU')
csvReader = csv.reader(f)
sheet = []
for row in csvReader:
    sheet.append(row)
data = [i for i in sheet[1:] if float(i.count('0'))/len(i) <= 0.1 and i[0]!='miRNA_ID']
data_1 = [[0 for i in range(len(data[0]))]for i in range(len(data))]
for i in range(len(data)):
    data[i][1:] = map(float,data[i][1:])
    for j in range(len(data[0])):
        if j == 0:
            data_1[i][j] = data[i][j]
        elif data[i][j] != 0:
            data_1[i][j] = math.log(data[i][j]/sum(data[i][1:])*1000000,2)
        else:
            data_1[i][j] = 'NaN'
#添加表头 
data_1.insert(0,sheet[0])
final_data = open('D:/clinical/TCGA-BLCA/miRNA.csv','wb')
csvWriter = csv.writer(final_data)
for row in data_1:
    csvWriter.writerow(row)
final_data.close()
'''


#def testResult()
def analyzeResult(label,model,DataVecs):
	proba = model.predict_proba(DataVecs)[:,1]
	predict = proba >= 0.5
	predict = predict.reshape((len(label),1))
	label = label.reshape((len(label),1))
	proba = proba.reshape((len(label),1))
	print ("\nAccuracy: %f %%" % (100. * sum(label == predict) / len(label)))
	print("Precision:%f\nRecall:%f\nF1score:%f\n" %(precision_score(label,predict),
	recall_score(label,predict),
	f1_score(label,predict)))
	try:
		result_auc = proba
		print ("\nRoc:%f\nAUPR:%f\n" % (roc_auc_score(label,result_auc),
			average_precision_score(label,result_auc)))
	except:
		print "ROC unavailable"
        
def analyzeResult_1(label,model,DataVecs,Result,actual,prediction):
    proba = model.predict_proba(DataVecs)[:,1]
    predict = proba >= 0.5
    np.set_printoptions(threshold=20)
    #print predict
    #print label
    #print proba
    #print roc_auc_score(label,proba)
    predict = predict.reshape((len(label),1))
    label = label.reshape((len(label),1))
    proba = proba.reshape((len(label),1))
    try:
        result_auc = proba
        Result.append(roc_auc_score(label,result_auc))
        actual.append(sum(label.tolist(),[]))
        prediction.append(sum(result_auc.tolist(),[]))
    except:
        print "ROC unavailable"
	
		
def balance_data(dataDataVecs,label):
	posdatavecs = dataDataVecs[label == 1]
	negdatavecs = dataDataVecs[label == 0]
	posdatavecs = posdatavecs[np.random.permutation(xrange(len(posdatavecs)))[0:len(negdatavecs)]]
	dataDataVecs = np.vstack((posdatavecs,negdatavecs))
	label = np.hstack([np.ones(len(negdatavecs)),np.zeros(len(negdatavecs))])
	return dataDataVecs,label

#将数据分成N组，准备做N折的cross-validation
def Group(X_pos,y_pos,X_neg,y_neg,N):
    X_set = []
    y_set = []
    X_1 = X_pos
    y_1 = y_pos
    X_2 = X_neg
    y_2 = y_neg
    for k in range(N):
        if k<N-1:
            X_train_pos,X_test_pos,y_train_pos,y_test_pos = train_test_split(X_1,y_1,test_size = float(1)/(N-k),random_state = None)
            X_train_neg,X_test_neg,y_train_neg,y_test_neg = train_test_split(X_2,y_2,test_size = float(1)/(N-k),random_state = None)
            X_set.append(np.append(X_test_pos,X_test_neg,axis=0))
            y_set.append(np.append(y_test_pos,y_test_neg,axis=0))
            X_1 = X_train_pos
            y_1 = y_train_pos
            X_2 = X_train_neg
            y_2 = y_train_neg
        else:
            X_set.append(np.append(X_train_pos,X_train_neg,axis=0))
            y_set.append(np.append(y_train_pos,y_train_neg,axis=0))
 
    return X_set,y_set

def feature_selection(X_train,y_train):
    #s = SelectFromModel(RandomForestRegressor())
    #s = SelectFromModel(ElasticNet(alpha = 1,l1_ratio = 0.5))
    s = SelectFromModel(LogisticRegression(penalty="l1", C=0.5))
    s.fit(X_train, y_train)
    s1 = s.get_support(indices=False)

    Result = [i for i in range(len(s1)) if s1[i]]
    Feature = [index[i] for i in Result]
    f1 = open('D:/clinical/TCGA-'+Cancer+'/TCGA-'+Cancer+'-'+molecule+'/'+drug+'-AUC/'+molecule+'_features_CN.csv','a')
    csvWriter = csv.writer(f1)
    csvWriter.writerow(Feature)
    f1.close()
    return Result
#五折的cross-validation
def cross_validation(X,y,n):
    Result = []
    actual = []
    prediction = []
    X_test = X[n]
    y_test = y[n]
    
    X_train = np.delete(X,n,axis = 0)
    y_train = np.delete(y,n,axis = 0)
    X_train_1 = X_train[0]
    y_train_1 = y_train[0]
    for j in range(1,3):
        X_train_1 = np.append(X_train_1,X_train[j],axis = 0)
        y_train_1 = np.append(y_train_1,y_train[j],axis = 0)
    features = feature_selection(X_train_1,y_train_1)

    X_train_1 = np.array([X_train_1.T[i] for i in features]).T
    X_test = np.array([X_test.T[i] for i in features]).T
    '''
    for i in range(len(X_train_1)):
        print X_train_1[i]
    print '\n'
    for i in range(len(X_test)):
        print X_test[i]
    '''
    scaler = preprocessing.StandardScaler().fit(X_train_1)
    scaler.transform(X_train_1)
    scaler.transform(X_test)
    for classifier in classifiers:
        classifier.fit(X_train_1, y_train_1)
        analyzeResult_1(y_test,classifier,X_test,Result,actual,prediction)
        
    return Result,actual,prediction
    



class LR(LogisticRegression):
    def __init__(self, threshold=0.01, dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None, solver='liblinear', max_iter=100,
                 multi_class='ovr', verbose=0, warm_start=False, n_jobs=1):

        #权值相近的阈值
        self.threshold = threshold
        LogisticRegression.__init__(self, penalty='l1', dual=dual, tol=tol, C=C,
                 fit_intercept=fit_intercept, intercept_scaling=intercept_scaling, class_weight=class_weight,
                 random_state=random_state, solver=solver, max_iter=max_iter,
                 multi_class=multi_class, verbose=verbose, warm_start=warm_start, n_jobs=n_jobs)
        #使用同样的参数创建L2逻辑回归
        self.l2 = LogisticRegression(penalty='l2', dual=dual, tol=tol, C=C, fit_intercept=fit_intercept, intercept_scaling=intercept_scaling, class_weight = class_weight, random_state=random_state, solver=solver, max_iter=max_iter, multi_class=multi_class, verbose=verbose, warm_start=warm_start, n_jobs=n_jobs)

    def fit(self, X, y, sample_weight=None):
        #训练L1逻辑回归
        super(LR, self).fit(X, y, sample_weight=sample_weight)
        self.coef_old_ = self.coef_.copy()
        #训练L2逻辑回归
        self.l2.fit(X, y, sample_weight=sample_weight)

        cntOfRow, cntOfCol = self.coef_.shape
        #权值系数矩阵的行数对应目标值的种类数目
        for i in range(cntOfRow):
            for j in range(cntOfCol):
                coef = self.coef_[i][j]
                #L1逻辑回归的权值系数不为0
                if coef != 0:
                    idx = [j]
                    #对应在L2逻辑回归中的权值系数
                    coef1 = self.l2.coef_[i][j]
                    for k in range(cntOfCol):
                        coef2 = self.l2.coef_[i][k]
                        #在L2逻辑回归中，权值系数之差小于设定的阈值，且在L1中对应的权值为0
                        if abs(coef1-coef2) < self.threshold and j != k and self.coef_[i][k] == 0:
                            idx.append(k)
                    #计算这一类特征的权值系数均值
                    mean = coef / len(idx)
                    self.coef_[i][idx] = mean
        return self
#调用的分类算法名称
names = ["DummyClassifier","Linear SVM", "RBF SVM",
		 "Decision Tree", "Random Forest",
		 "Naive Bayes"]
  
#实例化所用的分类算法，下面的参数可以根据官网API自行调整
classifiers = [
	DummyClassifier(strategy = "uniform"),
	SVC(kernel="linear", C=0.025,probability = True),
	SVC( C=1,probability = True),
	DecisionTreeClassifier(max_depth=None),
	RandomForestClassifier(max_depth=None, n_estimators=100, max_features='auto'),
	GaussianNB()]


#从TCGA数据中获取X和y
Cancer = 'CESC'
drug = 'Cisplatin'
molecule = 'mRNA'
gm = pd.read_table('D:/clinical/TCGA-'+Cancer+'/TCGA-'+Cancer+'-'+molecule+'/'+molecule+'.csv',sep = ",")
cd = pd.read_table('D:/clinical/TCGA-'+Cancer+'/TCGA-'+Cancer+'-CLINIC/TCGA-'+Cancer+'-'+drug+'.csv',sep= ",") 

#筛选有药物响应的TCGA编号的miRNA数据
#print gm.columns#列索引
#print gm.index#行索引
ID1 = gm.columns.values

ID2 = cd.SampleID.values

ID = list(set(ID1).intersection(set(ID2)))
print len(ID)
gm = gm.loc[:,ID]

#设定X就是基因表达量矩阵
X = np.array(gm).T

#处理缺省值
X = Imputer().fit_transform(X)
print len(X)
#查重及去重
data = DataFrame(cd)

'''
IsDuplicated = data.duplicated() 
for i in range(len(IsDuplicated)):
    print IsDuplicated[i]
'''
data = data.drop_duplicates()  

#根据sampleID索引对应的标签
data.index = data['SampleID']
y = (data['Response'][ID] == "Sensitive")

#L1正则化
#X = SelectFromModel(LogisticRegression(penalty="l1", C=0.3)).fit_transform(X, y)

#X = SelectFromModel(LR(threshold=0.4, C=0.5)).fit_transform(X, y)
#Feature Pre-selection,做t-test，筛选去p-value>0.01的
y1 = [i for i in range(len(y)) if y[i]]
y2 = [i for i in range(len(y)) if  not y[i]]

index = []
for k in range(len(X.T)):
    x1 = [X.T[k][i] for i in y1]
    x2 = [X.T[k][i] for i in y2]
    p_value_1 = levene(x1,x2)#检验是否有方差齐性
    if p_value_1 >0.05:
        p_value = ttest_ind(x1,x2,equal_var=True).pvalue
        if p_value <= 0.01:
            index.append(k)
    else:
        p_value = ttest_ind(x1,x2,equal_var=False).pvalue
        if p_value <= 0.01:
            index.append(k)
X = np.array([X.T[i] for i in index]).T
print index
print len(index)
#正负样本数量不一样，这边做balance问题，所以平衡正负样本数量

X,y = balance_data(X,y)
y = pd.Series(y)
'''
#对于每个数据集中的特征向量X和预测对象y
# 归一化特征并划分训练测试集
#X = StandardScaler().fit_transform(X)
'''
#正类样本所在位置
pos = [i for i in range(len(y)) if y[i]]
#负类样本所在位置
neg = [i for i in range(len(y)) if not y[i]] 

X_neg = [X[i] for i in neg]
X_pos = [X[i] for i in pos]
pos_ID = [y.index[i] for i in pos]
neg_ID = [y.index[i] for i in neg]
pos_bool = [y[i] for i in pos]
neg_bool = [y[i] for i in neg]
X_pos = np.array(X_pos)
X_neg = np.array(X_neg)
y_pos = pd.Series(pos_bool,index=pos_ID)
y_neg = pd.Series(neg_bool,index=neg_ID)

DummyClassifier = []
Linear_SVM = []
RBF_SVM = []
Decision_Tree = []
Random_Forest = []
Naive_Bayes = []
f = open('D:/clinical/TCGA-'+Cancer+'/TCGA-'+Cancer+'-'+molecule+'/'+drug+'-AUC/'+molecule+'-AUC-'+drug+'_CN.csv',"wb")
csvWriter = csv.writer(f)
csvWriter.writerow(['DummmyClassifer','Linear SVM','RBF SVM','Decision Tree','Random Forest','Naive Bayes'])
DummyClassifier_AC = []
Linear_SVM_AC = []
RBF_SVM_AC = []
Decision_Tree_AC = []
Random_Forest_AC = []
Naive_Bayes_AC = []
DummyClassifier_PRE = []
Linear_SVM_PRE = []
RBF_SVM_PRE = []
Decision_Tree_PRE = []
Random_Forest_PRE = []
Naive_Bayes_PRE = []
#for i in range(20):    
#五折cross-validation
#X_set,y_set=Group(X_pos,y_pos,X_neg,y_neg)
    
#c = np.append(X_set[0],X_set[1],axis = 0)
#a = np.concatenate([X_set[0],X_set[1]],axis = 0)
#b = np.delete(X_set,0,axis = 0)
for i in range(20):
    #5折cross-validation
    X_set,y_set=Group(X_pos,y_pos,X_neg,y_neg,5)
    for j in range(5):
        row = []
        Result_1,actual1,prediction1 = cross_validation(X_set,y_set,j)
        DummyClassifier.append(Result_1[0])
        Linear_SVM.append(Result_1[1])
        RBF_SVM.append(Result_1[2])
        Decision_Tree.append(Result_1[3])
        Random_Forest.append(Result_1[4])
        Naive_Bayes.append(Result_1[5])
        row.append(Result_1[0])
        row.append(Result_1[1])
        row.append(Result_1[2])
        row.append(Result_1[3])
        row.append(Result_1[4])
        row.append(Result_1[5])
        DummyClassifier_AC.extend(actual1[0])
        Linear_SVM_AC.extend(actual1[1])
        RBF_SVM_AC.extend(actual1[2])
        Decision_Tree_AC.extend(actual1[3])
        Random_Forest_AC.extend(actual1[4])
        Naive_Bayes_AC.extend(actual1[5])
        DummyClassifier_PRE.extend(prediction1[0])
        Linear_SVM_PRE.extend(prediction1[1])
        RBF_SVM_PRE.extend(prediction1[2])
        Decision_Tree_PRE.extend(prediction1[3])
        Random_Forest_PRE.extend(prediction1[4])
        Naive_Bayes_PRE.extend(prediction1[5])
        csvWriter.writerow(row)
f.close()
print np.array(DummyClassifier).mean()
print np.array(Linear_SVM).mean()
print np.array(RBF_SVM).mean()
print np.array(Decision_Tree).mean()
print np.array(Random_Forest).mean()
print np.array(Naive_Bayes).mean()
Classifier_AC = [DummyClassifier_AC,Linear_SVM_AC,RBF_SVM_AC,Decision_Tree_AC,Random_Forest_AC,Naive_Bayes_AC]
Classifier_PRE = [DummyClassifier_PRE,Linear_SVM_PRE,RBF_SVM_PRE,Decision_Tree_PRE,Random_Forest_PRE,Naive_Bayes_PRE]
#画出ROC曲线图
for name,AC,PRE in zip(names,Classifier_AC,Classifier_PRE):
    false_positive_rate, true_positive_rate, thresholds = roc_curve(AC, PRE)
    #plt.title(name)
    plt.plot(false_positive_rate, true_positive_rate, lw=1,
             label=name+'(AUC = %f)'%roc_auc_score(AC,PRE))
plt.legend(loc='lower right')
plt.plot([0,1],[0,1],'--')
plt.xlim([-0.1,1.2])
plt.ylim([-0.1,1.2])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.title('miRNA-Cisplatin-pan-cancer-random_AUC')
#plt.savefig('D:/clinical/TCGA-'+Cancer+'/TCGA-'+Cancer+'-'+molecule+'/'+drug+'-AUC/AUC_new.png')
plt.show()

#画出boxplot图
data = pd.DataFrame({'DymmyClassifier':DummyClassifier,
                     'Linear SVM':Linear_SVM,
                     'RBF SVM':RBF_SVM,
                     'Decision Tree':Decision_Tree,
                     'Random Forest':Random_Forest,
                     'Naive Bayes':Naive_Bayes})
data.boxplot()
plt.ylabel('AUC')
plt.xlabel('Classifier')
#plt.savefig('D:/clinical/TCGA-'+Cancer+'/TCGA-'+Cancer+'-'+molecule+'/'+drug+'-AUC/boxplot.png')
plt.show()

'''
#使用不同分类器进行训练/测试
for name, classifier in zip(names, classifiers):
	#使用训练集训练分类器
	print name
	classifier.fit(X_train, y_train)
	analyzeResult(y_test,classifier,X_test)
'''

    