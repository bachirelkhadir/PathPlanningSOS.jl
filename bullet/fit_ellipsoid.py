import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import random_projection
plt.style.use('dark_background')

COLLISION_MAP_FILE = "collsion-map-results-2021-01-08_21-46-27/collision_map.pd"
OUTPUT_FILE = os.path.join(os.path.dirname(COLLISION_MAP_FILE),
                           "collision_ellipsoid.pd")

OUTPUT_FILE
# "collsion-map-results-2021-01-06_22-34-12/collision_map.pd")
collision_map = pd.read_csv(COLLISION_MAP_FILE)
collision_map = collision_map[[f"j{i}" for i in range(7)] + ["collide"]]
eps =1e-2

training_idx = np.random.randint(len(collision_map), size=100000)
training_idx = list(range(len(collision_map)))
data = collision_map.iloc[training_idx, :7].values
y = collision_map.collide.values[training_idx]

C = np.cov(data[y].T)
Cinv = np.linalg.inv(C + eps*np.eye(7))
mu = np.mean(data[y].T, axis=1).round(2)

pd.DataFrame(np.concatenate([Cinv, [mu]]).T,
             columns=[ *[f"Q{i}" for i in range(7)], "mu" ])\
             .to_csv(OUTPUT_FILE)
##############################
# End
##############################


levelset = lambda u: (u - mu).T @ Cinv @ (u - mu).T
levelset_data = pd.DataFrame({"collide": y, "uQu": map(levelset, data)}).astype(float)
sns.set_theme()
sns_plot = sns.displot(levelset_data, x="uQu", hue="collide", kind="kde", )
sns_plot.set(xlim=(0, 50))




sns.heatmap(np.log(C).round(5), annot=True)

# dim reduction
#

# 2d
num_features = 2
scaler = StandardScaler().fit(data)
projection = PCA(num_features)  # project from 64 to 2 dimensions
# projection = random_projection.GaussianRandomProjection(num_features)
projected = projection.fit_transform(scaler.transform(data))

plt.scatter(projected[np.where(y), 0], projected[np.where(y), 1], c="yellow", s=2,  alpha=0.5)
plt.scatter(projected[np.where(~y), 0], projected[np.where(~y), 1], c="magenta", s=1,  alpha=0.5)
plt.axis("equal")


# 3d

num_features = 3
num_features
scaler = StandardScaler().fit(data)
projection = PCA(num_features).fit(scaler.transform(data))  # project from 64 to 2 dimensions
# projection = random_projection.GaussianRandomProjection(num_features)
projected = projection.transform(scaler.transform(collision_map.iloc[np.where(collision_map.collide)[0], :7]))

scatter_data = pd.DataFrame({"x": projected[:, 0],
                             "y": projected[:, 1],
                             "z": projected[:, 2],
                             "color": y,
                             "size": 1,
                             })
import plotly.express as px
fig = px.scatter_3d(scatter_data, x="x", y="y", z="z", color="color", size="size")
fig.write_html('collision_map_plotly.html', auto_open=False)



# nearest neighbour
from sklearn import neighbors
clf = neighbors.KNeighborsClassifier(1,)
clf.fit(data, y)

NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(data)
#


XX = collision_map.iloc[:, :7].values
X = XX[collision_map['collide']]



mu = np.mean(X, axis=0)
Sigma = np.linalg.inv(np.cov(X.T) + eps * np.eye(len(mu)))
lambda_, v = np.linalg.eig(Sigma)


shift = lambda u: u - mu
quad = lambda u: u.T @ Sigma @ u
f = lambda u: quad(shift(u))
XX.shape

XX_mu =XX[np.random.randint(len(XX), size=1000), :]
np.mean(list(map(f, XX_mu)))
X_mu =XX[np.random.randint(len(XX), size=1000), :]
np.mean(list(map(f, X_mu)))




X_sub = XX[np.random.randint(len(XX), size=1000), :]
pca = PCA(n_components=2)
pca.fit(X_sub)
X_pca = pca.transform(X_sub)
print("original shape:   ", X_sub.shape)
print("transformed shape:", X_pca.shape)
def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)

# plot data
plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.2)
for length, vector in zip(pca.explained_variance_, pca.components_):
    v = vector * 30 * np.sqrt(length)
    print(v)
    draw_vector(pca.mean_[:2], pca.mean_[:2] + v[:2])
plt.axis('equal');


plt.scatter(*XX_mu.T)
plt.scatter(*X_mu.T, color="r")


from matplotlib.patches import Ellipse

plt.figure()
ax = plt.gca()

ellipse = Ellipse(xy=(157.18, 68.4705), width=0.036, height=0.012,
                        edgecolor='r', fc='None', lw=2)
ax.add_patch(ellipse)




(XX_mu @ Sigma).shape



np.mean(list(map(f, X)))
np.mean(list(map(f, XX)))
