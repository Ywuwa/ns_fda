import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    # Trying to find module in the parent package
    from . import visualize
    print(True)
    
except ImportError:
    print('Relative import failed')

try:
    # Trying to find module on sys.path
    import visualize
    print(True)
except ModuleNotFoundError:
    print('Absolute import failed')

def get_frame(pfile):
  values = visualize.read_cube_txt(pfile)  # shape (nz, ny, nx)
  frame = visualize.visualize_nodes(values)
  return frame

num_frames = 300
fps = 12*2
fig, ax = plt.subplots()
ax.set_axis_off()  # отключим оси, чтобы было похоже на видео

# Создаем объект для записи видео через matplotlib
Writer = animation.FFMpegWriter
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)

with writer.saving(fig, "pressure_video.mp4", dpi=300):
  for i in range(num_frames):
    pfile = "p/" + f"{i}" + ".txt"
    frame = get_frame(pfile)
    # Отобразим кадр в `imshow`
    ax.imshow(frame)
    # Добавляем следующий кадр
    writer.grab_frame()
    # Очистим фигуру для следующего кадра
    ax.clear()
    ax.set_axis_off()
    if (i % 5 == 0):
      print(f"{i} of {num_frames}")

print("Видео сохранено как pressure_video.mp4")