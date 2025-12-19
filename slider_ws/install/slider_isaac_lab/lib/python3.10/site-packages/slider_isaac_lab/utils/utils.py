import matplotlib.pyplot as plt 
import torch 



def showImage(ax, image_torch : torch.Tensor): 
    image_numpy = image_torch.to("cpu").numpy()
    ax.clear()
    ax.imshow(image_numpy)
    ax.axis("off")
    plt.draw()
    plt.pause(0.1)