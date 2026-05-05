import cv2
import os
print(os.getcwd())

pasta_imagens = "/home/jvassalo/LBM-CUDA/animation/videos"
saida_video = "video.mp4"
fps = 30

# lista as imagens
arquivos = sorted(os.listdir(pasta_imagens))

# lê a primeira imagem para pegar o tamanho
frame = cv2.imread(os.path.join(pasta_imagens, arquivos[0]))
altura, largura, _ = frame.shape

# cria o vídeo
video = cv2.VideoWriter(
    saida_video,
    cv2.VideoWriter_fourcc(*'mp4v'),
    fps,
    (largura, altura)
)

# adiciona cada imagem
for nome in arquivos:
    caminho = os.path.join(pasta_imagens, nome)
    img = cv2.imread(caminho)
    video.write(img)

video.release()