from paraview.simple import *
import xml.etree.ElementTree as ET

tree = ET.parse("temp.xmf")
root = tree.getroot()

files = []

for inc in root.iter():
    if "include" in inc.tag:  # catches xi:include
        href = inc.attrib.get("href")
        if href:
            files.append(href)


def frender(name,view2,data,cloud,rain):
    
    # ============================================================
    # Data loading 
    # ============================================================
    data.FileNames = [name]
    data.UpdatePipeline()

    #apply thresholds
    cloud.UpdatePipeline()
    rain.UpdatePipeline()
    
    # ============================================================
    # Not a standard ParaView pipeline, but by recreating the display object 
    # quality of the volume rendering improves significantly.
    # ============================================================
    
    disp2 = Show(cloud, view2)
    disp2.Representation = 'Volume'
    disp2.ColorArrayName = ['CELLS', 'rc']
    disp2.InterpolateScalarsBeforeMapping = True

    # ============================================================
    # color (white cloud)
    # ============================================================
    lut = GetColorTransferFunction('rc')
    lut.RescaleTransferFunction(0.0, 0.004)
    lut.RGBPoints = [
        1e-5, 0.65, 0.70, 0.75,  # Cloud edge
        1e-4, 0.85, 0.88, 0.90,  # Cloud body
        1e-3, 1.00, 1.00, 1.00   # Thick convective core
    ]

    disp2.LookupTable = lut 

    #============================================================
    # opacity
    # ============================================================
    otf = GetOpacityTransferFunction('rc')
    otf.RescaleTransferFunction(0.0, 0.004)

    otf.Points = [
        1e-5, 0.10, 0.5, 0.0,
        1e-4, 0.60, 0.5, 0.0,
        1e-3, 0.80, 0.5, 0.0,
    ]
    disp2.OpacityTransferFunction = otf

    #============================================================
    # rendering tuning
    # ============================================================
    disp2.Ambient = 0.3
    disp2.Diffuse = 0.8
    disp2.Specular = 0.01
    disp2.SpecularPower = 10

    disp2.LookupTable = lut

    # Display rain as volume
    rain_disp = Show(rain, view2)
    rain_disp.Representation = 'Volume'

    ColorBy(rain_disp, ('CELLS', 'rr'))

    rain_lut = GetColorTransferFunction('rr')
    rain_lut.RGBPoints = [
    1e-7, 0.85, 0.90, 0.95,
    1e-5, 0.60, 0.75, 0.90,
    1e-4, 0.35, 0.55, 0.75,
    1e-3, 0.10, 0.20, 0.45,
    ]
    rain_disp.LookupTable = rain_lut
    print("Rendering frame: ", name[8:18])

    Render()
    SaveScreenshot("frames/frame"+name[8:18]+".png", view2)


#view setup
view2 = CreateView('RenderView')
view2.ViewSize = [600, 600]
view2.UseColorPaletteForBackground = 0
view2.BackgroundColorMode = "Gradient"
view2.Background  = [0.62, 0.70, 0.82]
view2.Background2 = [0.28, 0.36, 0.52]
view2.ResetCamera()
#camera setup
camera = view2.GetActiveCamera()
camera.SetFocalPoint(
    6020.906085968018,
    6020.90625,
    2446.76806640625
)

camera.SetViewUp(0.0, 0.0, 1.0)
camera.OrthogonalizeViewUp()
#camera.Azimuth(30)
#camera tilt
camera.Elevation(0)
#reader setup
data = XDMFReader(FileNames=files[0])
data.UpdatePipeline()
#threshold for rain, to avoid rendering zero values
cloud = Threshold(Input=data)
cloud.Scalars = ['CELLS', 'rc']
cloud.ThresholdMethod = 'Between'
cloud.LowerThreshold = 1e-9
cloud.UpperThreshold = 1
cloud.UpdatePipeline()
#threshold for rain, to avoid rendering small values
rain = Threshold(Input=data)
rain.Scalars = ['CELLS', 'rr']
rain.ThresholdMethod = 'Between'
rain.LowerThreshold = 1e-7
rain.UpperThreshold = 5e-3
rain.UpdatePipeline()
#Global Title, can be easily change every frame by modifying textproperty
title2 = Text()
title2.Text = "Cloud animation sample"
t2 = Show(title2, view2)
t2.WindowLocation = 'Upper Center'
t2.FontSize = 18
t2.Color = [0, 0, 0]
t2 = Show(title2, view2)

#outline = Outline(Input=data)
#Show(outline,view2)

#bypassing empty first frame
for fn in files[1:]:
    #camera rotation
    camera.Azimuth(0.4)
    
    frender(fn, view2,data,cloud,rain)

# after rendering use program in frames folder ffmpeg -framerate 24 -pattern_type glob        -i 'frame0000*.png'        -c:v libx264 -pix_fmt yuv420p  