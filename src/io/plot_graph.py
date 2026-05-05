import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("tke.csv")

plt.plot(data["time"], data["tke"])
plt.savefig("tke.png", dpi=300)
plt.clf()


data = pd.read_csv("vel.csv")

plt.plot(data["ux"], data["x"])
plt.xlabel("$u_x$")
plt.ylabel("$y/L$")
plt.grid()
plt.xlim(-1,1)
plt.ylim(0,1)
plt.axvline(x=0, color="black")
plt.savefig("ux.png", dpi=300)
plt.clf()


plt.plot(data["x"], data["uy"])
plt.ylabel("$u_y$")
plt.xlabel("$x/L$")
plt.xlim(0,1)
plt.ylim(-1,1)
plt.axhline(y=0, color="black")
plt.grid()
plt.savefig("uy.png", dpi=300)


