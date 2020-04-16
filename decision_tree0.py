import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns 
import matplotlib.pyplot as plt

from graphviz import Source
from pathlib import Path
import os

from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_graphviz
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix

def nad_tree_classifier(x, y): 
    tree_cf = DecisionTreeClassifier(max_depth=3)

    tree_cf.fit(x, y)

    return tree_cf

def save_fig(fig_id, tight_layout=True, fig_extension="png", resolution=300):
    path = os.path.join(IMAGES_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)

def main(): 

    PROJECT_ROOT_DIR = "."
    CHAPTER_ID = "decision_trees"
    IMAGES_PATH = os.path.join(PROJECT_ROOT_DIR, "images", CHAPTER_ID)
    os.makedirs(IMAGES_PATH, exist_ok=True)

    iris = load_iris()
    x = iris.data
    y = iris.target

    x_train, y_train, x_train, y_test = train_test_split(x,\
                                    y, test_size=0.2, random_state=42)

    tree_cf = nad_tree_classifier(x_train, y_train)

    y_pred = cross_val_predict(tree_cf, x_train, y_train, cv=10)

    print(y_pred)

    conf_matrix = confusion_matrix(y_train, y_pred)
    print(conf_matrix)

    # export_graphviz(
    #     tree_cf,
    #     out_file=os.path.join(IMAGES_PATH, "iris_tree.png"),
    #     feature_names=iris.feature_names[2:],
    #     class_names=iris.target_names,
    #     rounded=True,
    #     filled=True
    # )

    # Source.from_file(os.path.join(IMAGES_PATH, "iris_tree.png"))





if __name__ == "__main__":
    main()

