from thorcam.camera import ThorCam
from instrumental import list_instruments, instrument
import instrumental
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from instrumental.drivers.cameras import uc480

paramsets = list_instruments() ## camera found
print(paramsets)

paramsets = instrumental.list_instruments()
cammer = instrumental.instrument(paramsets[0])

fig1 = Figure(figsize=(6,6))
ax1 = fig1.add_subplot(111)
plt.figure()
framer= cammer.grab_image(timeout='1s',copy=True,n_frames=1,exposure_time='5ms',cx=640,
                                      left=10,cy=600,top=300)
plt.pcolormesh(framer)


# init camera
instruments = uc480.list_instruments()
cam = uc480.UC480_Camera(instruments[0])
# params
ls_frames = []
print(cam.DEFAULT_KWDS)
cam.start_live_video(framerate = "30Hz")
print (cam.is_open, cam.framerate)
#cam._set_exposure("5ms")
#cam._set_AOI(0, 0, 180, 120)

while (i < total_frames):
    flag = cam.wait_for_frame(timeout = "0 ms")
    if flag:
        i = i + 1
        frame = cam.latest_frame(copy = False)
        frame = frame.astype('uint8')
        ls_frames.append(frame)


class MyThorCam(ThorCam):
    def received_camera_response(self, msg, value):
        super(MyThorCam, self).received_camera_response(msg, value)
        if msg == 'image':
            return
        print('Received "{}" with value "{}"'.format(msg, value))
    def got_image(self, image, count, queued_count, t):
        print('Received image "{}" with time "{}" and counts "{}", "{}"'
              .format(image, t, count, queued_count))

cam = MyThorCam()
cam.start_cam_process()
cam.refresh_cameras()
