# 00. Module importing ====================================================
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from scipy import stats
import multiprocessing as mp
import os
# ========================================================================

# 01. Training and test data importing ===================================
Train = pd.read_csv("CD238.trn.CD44", sep=" ").T
#print(Train)
Train["Group"] = Train.index.to_series().str.contains('KDC').map(lambda x: 0 if x else 1)
#Train.columns = Train.iloc[0, :]
#Train = Train.iloc[2:, :]
Train.columns=list(Train.columns[0:-1])+["Group"]
print(Train.iloc[0:5,0:5])
Test = pd.read_csv("CD238.tst.CD44", sep=" ").T
Test["Group"] = Test.index.to_series().str.contains('KDC').map(lambda x: 0 if x else 1)
#Test.columns = Test.iloc[0, :]
#Test = Test.iloc[2:, :]
Test.columns=list(Test.columns[0:-1])+["Group"]
print(Test.iloc[0:5,0:5])
feature_cols = [col for col in Train.columns.intersection(Test.columns) if col != "Group"]
X = Train[feature_cols]
y = Train["Group"]
X_test = Test[feature_cols]
y_test = Test["Group"]
print("#===============Data importing process completed================#")

# ========================================================================

# 02. 데이터 스케일링 (SVM, LogisticRegression만 적용)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X_test_scaled = scaler.transform(X_test)

print("#===============Data scaling completed================#")
# ========================================================================
from sklearn.utils.class_weight import compute_class_weight
import numpy as np

classes = np.unique(y)
weights = compute_class_weight(class_weight='balanced', classes=classes, y=y)
class_weight_dict = dict(zip(classes, weights))


# 03. 모델 및 파라미터 정의
models = {
    'LogisticRegression': LogisticRegression(random_state=42, max_iter=2000),
    'RandomForest': RandomForestClassifier(random_state=42, n_jobs=-1),
    'SVM': SVC(random_state=42),  # probability 제거
    'GradientBoosting': GradientBoostingClassifier(random_state=42)
}

param_grids = {
      'LogisticRegression':[{
         'penalty':[None],      # penalty 완전 제거 [web:10][web:16]
         'solver':['lbfgs'],    # liblinear 불필요, lbfgs가 안정적 [cite:7]
         'max_iter':[1000]}],     # 수렴 1000

      'RandomForest': {         # ★ 가장 추천
         'n_estimators': [500],
         'max_features': ['sqrt'],
         'min_samples_leaf': [5, 10],
         'min_samples_split': [10],
         'min_impurity_decrease': [0.001]
      },
      'GradientBoosting': {
         'n_estimators': [100,200,300],
         'learning_rate': [0.01,0.1,0.2],
         'subsample': [0.8,1.0],
         'max_features': ['sqrt'],
         'min_samples_leaf': [1,2]
      },
      'SVM': {
        'C': [0.1],
        'kernel': ['linear']
      } }


#param_grids = {
#     'LogisticRegression':[{
#        'penalty':[None],      # penalty 완전 제거 [web:10][web:16]
#        'solver':['lbfgs'],    # liblinear 불필요, lbfgs가 안정적 [cite:7]
#        'max_iter':[1000]}],     # 수렴 1000
#    'RandomForest': {           # ★ 45개 SNP용으로 확장
#        'n_estimators': [200, 500, 800],
#        'max_features': [None,0.8],
#        'min_samples_leaf': [1, 3, 5, 10],
#        'min_samples_split': [2, 5, 10],
#        'max_depth': [10, 20, None],
#        'min_impurity_decrease': [0.0, 0.001]
#    },
#
#    'GradientBoosting': {       # ★ 확장
#        'n_estimators': [100, 200, 500],
#        'learning_rate': [0.01, 0.05, 0.1, 0.2],
#        'subsample': [0.7, 0.8, 0.9, 1.0],
#        'max_features': [None, 0.8],
#        'max_depth': [3, 5, 7],
#        'min_samples_leaf': [1, 2, 5]
#    },
#    'SVM': {
#        'C': [0.1],
#        'kernel': ['linear']
#    }
#}

#param_grids = {
#    'LogisticRegression':[{
#        'penalty':[None],      # penalty 완전 제거 [web:10][web:16]
#        'solver':['lbfgs'],    # liblinear 불필요, lbfgs가 안정적 [cite:7]
#        'max_iter':[1000]}],     # 수렴 1000
#    'RandomForest': {
#        'n_estimators': [100, 200, 300, 500],
#        'max_depth': [None, 10, 20],
#        'min_samples_split': [2, 5],
#        'min_samples_leaf': [1, 2],
#        'max_features': [None],
#        'bootstrap': [True],
#        'class_weight': [None, 'balanced'],  # dict 제거
#        'criterion': ['gini', 'entropy'],  # log_loss 제거 (버전에 따라 에러)
#        'max_leaf_nodes': [None, 500],
#        'min_impurity_decrease': [0.0, 0.001],
#        'ccp_alpha': [0.0, 0.001],
#        'max_samples': [None, 0.8]  # stochasticity 강화
#    },
#    'SVM': {
#        'C': [0.1, 1, 10],
#        'kernel': ['rbf', 'linear'],
#        'gamma': ['scale', 0.01, 0.1],  # poly 제거로 충돌 방지
#        'class_weight': [None, 'balanced'],  # dict 제거
#        'shrinking': [True],
#        'max_iter': [-1],  # 무제한
#        'tol': [1e-3]
#    },
#    'GradientBoosting': {
#        'n_estimators': [100, 200, 300],
#       'learning_rate': [0.01, 0.1, 0.2],
#       'max_depth': [3, 5, 7],
#       'subsample': [0.8, 1.0],
#       'min_samples_split': [2, 5],
#       'min_samples_leaf': [1, 2],
#       'max_features': [None]  # auto 제거 (deprecated)
#   }
#}


print("#===Model definition process is end===#")
# ========================================================================

# 04. 모델별 데이터 매핑
X_dict = {
    'LogisticRegression': X_scaled,
    'SVM': X_scaled,
    'RandomForest': X,
    'GradientBoosting': X
}
X_test_dict = {
    'LogisticRegression': X_test_scaled,
    'SVM': X_test_scaled,
    'RandomForest': X_test,
    'GradientBoosting': X_test
}

print("#===Model setting process is end===#")
# ========================================================================

# 05. 하이퍼파라미터 튜닝 및 학습 함수
from sklearn.model_selection import StratifiedKFold
kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

def train_and_evaluate_model(model, param_grid, X_train, y_train, X_test, y_test):
    if isinstance(model, (RandomForestClassifier, GradientBoostingClassifier)):
        search = RandomizedSearchCV(model, param_grid, n_iter=200, cv=kfold, scoring='roc_auc', random_state=42)
    else:
        search = GridSearchCV(model, param_grid, cv=kfold, scoring='roc_auc')
    search.fit(X_train, y_train)
    best_model = search.best_estimator_

    if hasattr(best_model, "predict_proba"):
        y_pred_proba = best_model.predict_proba(X_test)[:, 1]
    else:
        y_pred_proba = best_model.decision_function(X_test)

    thresholds = np.arange(0.01, 0.99, 0.01)
    best_f1 = 0
    best_threshold = 0.5
    for t in thresholds:
        preds = (y_pred_proba >= t).astype(int)
        score = f1_score(y_test, preds)
        if score > best_f1:
            best_f1 = score
            best_threshold = t

    y_pred = (y_pred_proba >= best_threshold).astype(int)
    auc_roc = roc_auc_score(y_test, y_pred_proba)
    f1 = f1_score(y_test, y_pred)
    n = len(y_test)
    se_auc = np.sqrt(auc_roc * (1 - auc_roc) / n)
    ci_auc_lower, ci_auc_upper = stats.norm.interval(0.95, loc=auc_roc, scale=se_auc)
    return auc_roc, f1, (ci_auc_lower, ci_auc_upper), search.best_params_, best_threshold, y_pred_proba, y_pred

print("#===Model tunning and hyperparameter setting process is end===#")
# ========================================================================

# 06. 병렬 처리 및 모델별 학습
def parallel_train_and_evaluate(args):
    model_name, model, param_grid, X_train, y_train, X_test, y_test = args
    auc_roc, f1, ci_auc, best_params, best_threshold, y_pred_proba, y_pred = train_and_evaluate_model(
        model, param_grid, X_train, y_train, X_test, y_test
    )
    return model_name, auc_roc, f1, ci_auc, best_params, best_threshold, y_pred_proba, y_pred

args_list = [
    (model_name, models[model_name], param_grids[model_name],
     X_dict[model_name], y, X_test_dict[model_name], y_test)
    for model_name in models.keys()
]

with mp.Pool(processes=4) as pool:
    parallel_results = pool.map(parallel_train_and_evaluate, args_list)

results = {}
sample_pred_dict = {}
for model_name, auc_roc, f1, ci_auc, best_params, best_threshold, y_pred_proba, y_pred in parallel_results:
    results[model_name] = {
        'auc_roc': auc_roc,
        'f1': f1,
        'ci_auc': ci_auc,
        'best_params': best_params,
        'threshold': best_threshold,
        'y_pred_proba': y_pred_proba,
        'y_pred': y_pred
    }
    sample_pred_dict[model_name] = {
        "prs": y_pred_proba,
        "pred_class": y_pred,
        "threshold": best_threshold
    }

print("#===Model setting process is end, hyperparameter tuning done===#")
# ========================================================================

# 07. 최적 파라미터로 모델 재정의 (SVM은 probability 파라미터 제거)
svm_params = results['SVM']['best_params'].copy()
svm_params.pop('probability', None)
models_final = {
    'LogisticRegression': LogisticRegression(**results['LogisticRegression']['best_params'], random_state=42),
    'RandomForest': RandomForestClassifier(**results['RandomForest']['best_params'], random_state=42),
    'SVM': SVC(**svm_params, random_state=42),
    'GradientBoosting': GradientBoostingClassifier(**results['GradientBoosting']['best_params'], random_state=42)
}

print("#===Model processing...===#")
# ========================================================================

# 08. 결과를 DataFrame으로 변환 및 파일로 저장 ============================
output_name = "output/CD44_260303.log_default_v5"
os.makedirs(os.path.dirname(output_name), exist_ok=True)
RSID_1 = [x for x in Train.columns if x.startswith("rs")]
RSID_1 = '+'.join(RSID_1)

results_df = pd.DataFrame({
    'combination': [RSID_1],
    **{f'{model_name}_auc_roc': [scores['auc_roc']] for model_name, scores in results.items()},
    **{f'{model_name}_f1': [scores['f1']] for model_name, scores in results.items()},
    **{f'{model_name}_ci_auc_lower': [scores['ci_auc'][0]] for model_name, scores in results.items()},
    **{f'{model_name}_ci_auc_upper': [scores['ci_auc'][1]] for model_name, scores in results.items()},
    **{f'{model_name}_best_params': [scores['best_params']] for model_name, scores in results.items()},
    **{f'{model_name}_best_threshold': [scores['threshold']] for model_name, scores in results.items()}
})

results_df['combo_length'] = [X.shape[1]]

output_file = output_name + "_all_models_metrics_with_hyperparams.260202.txt"
results_df.to_csv(output_file, index=False, sep="\t")

# 09. 샘플별 예측값 및 PRS 평균/합산 저장 ================================
sample_preds = pd.DataFrame({
    "Sample": X_test.index,
    "TrueClass": y_test.astype(int).values
})

# 각 모델별 PRS, 예측 클래스, threshold 추가
for model_name in models.keys():
    sample_preds[f"{model_name}_PRS"] = sample_pred_dict[model_name]["prs"]
    sample_preds[f"{model_name}_PredClass"] = sample_pred_dict[model_name]["pred_class"]
    sample_preds[f"{model_name}_Threshold"] = sample_pred_dict[model_name]["threshold"]

# PRS 합산 및 평균
prs_cols = [f"{model_name}_PRS" for model_name in models.keys()]
sample_preds["PRS_sum"] = sample_preds[prs_cols].sum(axis=1)
sample_preds["PRS_mean"] = sample_preds[prs_cols].mean(axis=1)

sample_preds.to_csv(output_name + "_sample_predictions.260202.txt", sep="\t", index=False)

print("모든 결과 저장 완료!")


# ★★★ 10. SHAP 분석 및 개별 모델 시각화 (모든 모델) ★★★
import shap
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

print("#===SHAP analysis for ALL models started===#")

# 커스텀 SHAP scatter plot 함수 정의
def plot_shap_scatter(model_name, shap_values, X_data, feature_cols, output_name, top_n=15):
    """개별 모델별 SHAP scatter + summary plot"""

    plt.figure(figsize=(20, 6))

    # 1) Top N SHAP Bar plot
    plt.subplot(1, 3, 1)
    shap_importance = np.abs(shap_values).mean(0)
    top_indices = np.argsort(shap_importance)[-top_n:]
    top_snps = [feature_cols[i] for i in top_indices]
    top_shap = shap_importance[top_indices]

    sns.barplot(x=top_shap, y=top_snps, palette='viridis')
    plt.title(f'{model_name}\nTop {top_n} SNPs (AUC: {results[model_name]["auc_roc"]:.3f})')
    plt.xlabel('Mean |SHAP value|')

    # 2) SHAP Summary (beeswarm)
    plt.subplot(1, 3, 2)
    shap.summary_plot(shap_values, X_data, feature_names=feature_cols, 
                      max_display=15, show=False, plot_type="dot")
    plt.title(f'{model_name} SHAP Summary')

    # 3) Top1 SNP SHAP Dependence scatter
    plt.subplot(1, 3, 3)
    top_snp_idx = top_indices[0]
    top_snp_name = feature_cols[top_snp_idx]
    shap.plots.scatter(shap_values[:, top_snp_idx], 
                      color=shap_values[:, top_indices[1]],  # 2nd top으로 색상
                      show=False, alpha=0.6)
    plt.xlabel(top_snp_name)
    plt.ylabel('SHAP value')
    plt.title(f'{model_name}: {top_snp_name} dependence')

    plt.tight_layout()
    plt.savefig(f'{output_name}_{model_name}_shap_complete.png', dpi=300, bbox_inches='tight')
    plt.close()

    return top_indices, shap_importance

# 각 모델별 SHAP 분석 실행
shap_results = {}
for model_name in models.keys():
    print(f"Analyzing {model_name}...")

    # 모델별 적절한 explainer 선택
    if model_name == 'SVM':
        explainer = shap.LinearExplainer(models_final[model_name], X_test_dict[model_name])
    elif model_name == 'LogisticRegression':
        explainer = shap.LinearExplainer(models_final[model_name], X_test_dict[model_name])
    else:  # RF, GBM (TreeExplainer)
        explainer = shap.TreeExplainer(models_final[model_name])

    shap_values = explainer.shap_values(X_test_dict[model_name])

    # Scatter plot 생성
    top_indices, shap_importance = plot_shap_scatter(
        model_name, shap_values, X_test_dict[model_name], 
        feature_cols, output_name
    )

    # Top 10 SNP 저장
    top_snps = [feature_cols[i] for i in top_indices[:10]]
    shap_results[model_name] = {
        'top10_snps': top_snps,
        'shap_importance': shap_importance,
        'shap_values': shap_values
    }

    print(f"{model_name} Top 3: {top_snps[:3]}")

# 4개 모델 비교: Top SNP 일관성 확인
print("\n=== MODEL COMPARISON ===")
comparison_df = pd.DataFrame({
    'Model': list(models.keys()),
    'Top1_SNP': [shap_results[m]['top10_snps'][0] for m in models.keys()],
    'AUC': [results[m]['auc_roc'] for m in models.keys()]
})
print(comparison_df)

# 모든 모델 공통 Top 10 SNPs 찾기
all_top_snps = []
for model_name in models.keys():
    all_top_snps.extend(shap_results[model_name]['top10_snps'][:5])
common_snps = pd.Series(all_top_snps).value_counts().head(10)
print(f"\n공통 Top 10 SNPs:\n{common_snps}")

# 공통 SNP Barplot
plt.figure(figsize=(12, 6))
common_snps.plot(kind='bar')
plt.title('Common Top SNPs Across All Models')
plt.ylabel('Appearance Count (out of 4 models)')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f'{output_name}_common_top_snps.png', dpi=300, bbox_inches='tight')
plt.close()

# SHAP 결과 저장
shap_summary = pd.DataFrame({
    'Model': list(models.keys()),
    'Top1_SNP': [shap_results[m]['top10_snps'][0] for m in models.keys()],
    'AUC': [results[m]['auc_roc'] for m in models.keys()],
    **{f'Top{k+1}_SNP': [shap_results[m]['top10_snps'][k] for m in models.keys()] 
       for k in range(5)}
})
shap_summary.to_csv(f'{output_name}_shap_model_comparison.txt', sep='\t', index=False)

print("#===============SHAP analysis for ALL models completed================#")


