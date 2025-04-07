import argparse
import sys
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def parse_args():
    parser = argparse.ArgumentParser(description="Load an XMF file and settings.")
    parser.add_argument("file_path", type=str, help="Path to the XMF file.")
    parser.add_argument(
        "--real",
        action="store_true",
        help="Changes the visualization to a more realistic one, with a background image and without the gridlines.",
    )
    parser.add_argument(
        "--lowerthreshold", type=float, help="Lower threshold value.", default=0.0001
    )
    parser.add_argument(
        "--boxcolor",
        type=float,
        nargs=3,
        metavar=("R", "G", "B"),
        default=[0.0, 0.0, 0.0],
        help="Box color (R G B) in range [0,1].",
    )
    parser.add_argument(
        "--boxsize",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=[1.0, 1.0, 0.6],
        help="Box size (X Y Z) in km.",
    )
    parser.add_argument(
        "--background",
        type=float,
        nargs=3,
        metavar=("R", "G", "B"),
        default=[0.67, 1.0, 1.0],
        help="Background color (R G B) in range [0,1].",
    )
    parser.add_argument(
        "--photopath", type=str, help="Path to the backgroung photo (.jpg)."
    )
    parser.add_argument(
        "--colormap1",
        type=str,
        default="Black-Body Radiation",
        help="Colormap for Calculator1.",
    )
    parser.add_argument(
        "--colormap2",
        type=str,
        default="Rainbow Uniform",
        help="Colormap for Calculator2.",
    )
    parser.add_argument(
        "--opacity1", type=float, default=1.0, help="Opacity for Calculator1."
    )
    parser.add_argument(
        "--opacity2", type=float, default=0.17, help="Opacity for Calculator2."
    )
    parser.add_argument(
        "--logscale1", action="store_true", help="Use log scale for Calculator1."
    )
    parser.add_argument(
        "--logscale2", action="store_true", help="Use log scale for Calculator2."
    )
    parser.add_argument(
        "--contour1",
        action="store_true",
        default=False,
        help="Use contour for Calculator1.",
    )
    parser.add_argument(
        "--contour2",
        action="store_true",
        default=False,
        help="Use contour for Calculator2.",
    )
    parser.add_argument(
        "--camera",
        type=str,
        choices=["default", "X+", "X-", "Y+", "Y-", "Z+", "Z-"],
        help="Camera placement option.",
        default="default",
    )
    parser.add_argument(
        "--middleframename",
        type=str,
        help="Name of the file with middle frame",
        default="anim_frame_middle_time.pdf",
    )
    parser.add_argument(
        "--animationname",
        type=str,
        help="Name of the file with animation",
        default="animationuwrealistic.ogv",
    )
    parser.add_argument(
        "--framerate", type=int, help="Number of frame rates.", default=5
    )
    parser.add_argument(
        "--axisinrealistic",
        action="store_true",
        help="Add axis to realistic visualization.",
    )
    return parser.parse_args(sys.argv[1:])


def load_xmf(file_path):
    tempxmf = XDMFReader(registrationName="temp.xmf", FileNames=[file_path])
    return tempxmf


def display_settings(red, green, blue, x=1.0, y=1.0, z=0.6):
    tempxmfDisplay = Show(tempxmf, renderView1, "StructuredGridRepresentation")
    tempxmfDisplay.Representation = "Outline"
    renderView1.Update()
    animationScene1.AnimationTime = 7600.0
    tempxmfDisplay.AmbientColor = [red, green, blue]
    tempxmfDisplay.DiffuseColor = [red, green, blue]
    tempxmfDisplay.Scale = [x, y, z]
    tempxmfDisplay.DataAxesGrid.Scale = [x, y, z]
    tempxmfDisplay.PolarAxes.Scale = [x, y, z]


def axes_settings():
    renderView1.AxesGrid.Visibility = 1
    renderView1.AxesGrid.XTitle = "X [km]"
    renderView1.AxesGrid.YTitle = "Y [km]"
    renderView1.AxesGrid.ZTitle = "Z [km]"
    renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.DataScale = [1000.0, 1000.0, 1000.0]
    renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ShowGrid = 1
    renderView1.AxesGrid.DataScale = [1000.0, 1000.0, 800.0]
    renderView1.ViewSize = [1000, 800]
    renderView1.AxesGrid.ZAxisUseCustomLabels = 1
    renderView1.AxesGrid.ZAxisLabels = [0.0, 2.0, 4.0, 6.0, 8.0]
    renderView1.AxesGrid.XAxisUseCustomLabels = 1
    renderView1.AxesGrid.XAxisLabels = [0.0, 2.0, 4.0, 6.0, 8, 0, 10.0, 12.0]
    renderView1.AxesGrid.YAxisUseCustomLabels = 1
    renderView1.AxesGrid.YAxisLabels = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
    renderView1.AxesGrid.CustomBounds = [0.0, 12000.0, 0.0, 12000.0, 0.0, 8000]
    renderView1.AxesGrid.UseCustomBounds = 1


def background(pathtophoto=None, red=0.67, green=1.0, blue=1.0, photobackgroung=False):
    renderView1.UseColorPaletteForBackground = 0
    renderView1.Background = [red, green, blue]
    if photobackgroung is True:
        texture = servermanager.rendering.ImageTexture()
        texture.FileName = pathtophoto
        GetActiveView().BackgroundColorMode = "Texture"
        GetActiveView().BackgroundTexture = texture
        Render()


def create_treshold(
    registrationName, input1, representation, lowertreshold, scalars=False
):
    threshold = Threshold(registrationName=registrationName, Input=input1)
    thresholdDisplay = Show(threshold, renderView1, representation)
    threshold.LowerThreshold = lowertreshold
    if scalars is True:
        threshold.Scalars = ["CELLS", "rr"]
    Hide(threshold, renderView1)
    return threshold


def calculator(
    registrationname,
    input,
    function,
    resultarrayname,
    colormap,
    scalarbartitle,
    opacity,
    x_move=0.0,
    y_move=0.0,
    z_move=7000.0,
    xscalarbar=None,
    yscalarbar=None,
    logscale=True,
    contour_value=False,
    move_data=False,
):
    calculator = Calculator(registrationName=registrationname, Input=input)
    calculatorDisplay = Show(calculator, renderView1, "UnstructuredGridRepresentation")
    calculatorDisplay.Representation = "Surface"
    calculatorDisplay.SetScalarBarVisibility(renderView1, True)
    calculator.Function = function
    calculatorDisplay.Opacity = opacity
    calculator.ResultArrayName = resultarrayname
    renderView1.Update()
    ColorBy(calculatorDisplay, ("CELLS", resultarrayname))
    calculatorDisplay.RescaleTransferFunctionToDataRange(True, False)
    rgkgLUT = GetColorTransferFunction(resultarrayname)
    rgkgLUT.ApplyPreset(colormap, True)
    rgkgLUT.MapControlPointsToLogSpace()
    if logscale is True:
        rgkgLUT.UseLogScale = 1
    if contour_value is False:
        rgkgLUTColorBar = GetScalarBar(rgkgLUT, renderView1)
        rgkgLUTColorBar.ComponentTitle = ""
        rgkgLUTColorBar.Title = scalarbartitle
        rgkgLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        rgkgLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        rgkgLUTColorBar.WindowLocation = "Any Location"
        rgkgLUTColorBar.Position = [xscalarbar, yscalarbar]
    if contour_value is True:
        ColorBy(calculatorDisplay, None)
        calculatorDisplay = GetDisplayProperties(calculator, view=renderView1)
        renderView1.OrientationAxesXVisibility = 0
        renderView1.OrientationAxesYVisibility = 0
        renderView1.OrientationAxesZVisibility = 0
        renderView1.OrientationAxesVisibility = 0
        renderView1.AxesGrid.Visibility = 0
    if move_data is True:
        ColorBy(calculatorDisplay, None)
        transform = Transform(registrationName="Transform", Input=calculator)
        transformDisplay = Show(
            transform, renderView1, "UnstructuredGridRepresentation"
        )
        transformDisplay.Representation = "Surface"
        transformDisplay.Opacity = opacity
        Hide(calculator, renderView1)
        transform.Transform.Translate = [x_move, y_move, z_move]
        renderView1.Update()


def hide_scalar_bar():
    rcLUT = GetColorTransferFunction("rc")
    HideScalarBarIfNotNeeded(rcLUT, renderView1)
    rrLUT = GetColorTransferFunction("rr")
    HideScalarBarIfNotNeeded(rrLUT, renderView1)


def time(view, position_x, position_y):
    time1 = AnnotateTimeFilter(guiName="AnnotateTimeFilter1", Format="Time:{time:f}s")
    timeDisplay = Show(time1, view)
    timeDisplay.WindowLocation = "Any Location"
    timeDisplay.Position = [position_x, position_y]
    timeDisplay.Color = [0.0, 0.0, 0.0]
    renderView1.Update()


def layout():
    layout1 = GetLayout()
    layout1.SetSize(1205, 749)


def current_camera_placement(choice):
    camera_settings = {
        "default": {
            "CameraPosition": [
                32031.133219522875 + 1000,
                -23188.99873874442 - 500,
                8278.119594860653 + 1000,
            ],
            "CameraFocalPoint": [
                6020.906250000001,
                6020.906250000002,
                6520.912597656248,
            ],
            "CameraViewUp": [0.0, 0.0, 1.0],
            "CameraParallelScale": 9136.058930334293,
        },
        "X+": {
            "CameraPosition": [-29279.811220785792, 6020.90625, 3312.5475585937497],
            "CameraViewUp": [0.0, 0.0, 1.0],
        },
        "X-": {
            "CameraPosition": [41321.62372078579, 6020.90625, 3312.5475585937497],
            "CameraViewUp": [0.0, 0.0, 1.0],
        },
        "Y+": {
            "CameraPosition": [6020.90625, -29279.811220785792, 3312.5475585937497],
            "CameraViewUp": [0.0, 1.0, 0.0],
        },
        "Y-": {
            "CameraPosition": [6020.90625, 41321.62372078579, 3312.5475585937497],
            "CameraViewUp": [0.0, 1.0, 0.0],
        },
        "Z+": {
            "CameraPosition": [6020.90625, 6020.90625, -31988.169912192043],
            "CameraViewUp": [1.0, 0.0, 0.0],
        },
        "Z-": {
            "CameraPosition": [6020.90625, 6020.90625, 38613.265029379545],
            "CameraViewUp": [1.0, 0.0, 0.0],
        },
        "Zoom": {
            "CameraPosition": [18000, -8000, 7000],
            "CameraFocalPoint": [6020.91, 6020.91, 6520.91],
            "CameraViewUp": [0.0, 0.0, 1.0],
            "CameraParallelScale": 1500.0,
        },
    }
    if choice not in camera_settings:
        raise ValueError("Invalid option, please enter a correct value.")
    renderView1.CameraPosition = camera_settings[choice].get(
        "CameraPosition", [6020.90625, 6020.90625, 3312.5475585937497]
    )
    renderView1.CameraFocalPoint = camera_settings["default"].get(
        "CameraFocalPoint", [6020.90625, 6020.90625, 3312.5475585937497]
    )
    renderView1.CameraViewUp = camera_settings[choice].get(
        "CameraViewUp", [0.0, 1.0, 0.0]
    )
    renderView1.CameraParallelScale = camera_settings["default"].get(
        "CameraParallelScale", 9136.49798722265
    )


def extract_bottom_xy(up_and_down=5000.0):
    renderView1.AxesGrid.Visibility = 1
    renderView1.AxesGrid.AxesToLabel = 2
    renderView1.AxesGrid.AxesToLabel = 10
    renderView1.AxesGrid.FacesToRender = 4
    renderView1.AxesGrid.XLabelFontSize = 20
    renderView1.AxesGrid.XTitleFontSize = 20
    renderView1.AxesGrid.XTitle = "X [km]"
    renderView1.AxesGrid.DataScale = [1000.0, 1000.0, 1000.0]
    renderView1.AxesGrid.XAxisUseCustomLabels = 1
    renderView1.AxesGrid.XAxisLabels = [0.0, 2.0, 4.0, 6.0, 8, 0, 10.0, 12.0]
    renderView1.AxesGrid.YLabelFontSize = 20
    renderView1.AxesGrid.YTitleFontSize = 20
    renderView1.AxesGrid.YTitle = "Y [km]"
    renderView1.AxesGrid.YAxisUseCustomLabels = 1
    renderView1.AxesGrid.YAxisLabels = [0.0, 2.0, 4.0, 6.0, 8, 0, 10.0, 12.0]
    renderView1.AxesGrid.CustomBounds = [0.0, 12000.0, 0.0, 12000.0, up_and_down, 8000]
    renderView1.AxesGrid.UseCustomBounds = 1


def save_middle_frame(filename):
    time_steps = tempxmf.TimestepValues
    middle_time = time_steps[len(time_steps) // 2]
    renderView1.ViewTime = middle_time
    for reader in (tempxmf,):
        reader.UpdatePipeline(middle_time)
        ExportView(
            filename=filename,
            view=renderView1,
            Rasterize3Dgeometry=False,
            GL2PSdepthsortmethod="BSP sorting (slow, best)",
        )
    RenderAllViews()


if __name__ == "__main__":
    args = parse_args()
    tempxmf = load_xmf(args.file_path)
    animationScene1 = GetAnimationScene()
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    renderView1 = GetActiveViewOrCreate("RenderView")
    materialLibrary1 = GetMaterialLibrary()
    if args.real is False:
        display_settings(*args.boxcolor, *args.boxsize)
        axes_settings()
        background()
        threshold1 = create_treshold(
            "Threshold1", tempxmf, "UnstructuredGridRepresentation", args.lowerthreshold
        )
        threshold2 = create_treshold(
            "Threshold2",
            tempxmf,
            "UnstructuredGridRepresentation",
            args.lowerthreshold,
            scalars=True,
        )
        calculator(
            registrationname="Calculator1",
            input=threshold1,
            function="rc*1000",
            resultarrayname="rc [g/kg]",
            colormap=args.colormap1,
            scalarbartitle="cloud water mixing ratio [g/kg]",
            opacity=args.opacity1,
            xscalarbar=0.9,
            yscalarbar=0.6,
            logscale=args.logscale1,
            contour_value=args.contour1,
        )
        calculator(
            registrationname="Calculator2",
            input=threshold2,
            function="rr*1000",
            resultarrayname="rr [g/kg]",
            colormap=args.colormap2,
            scalarbartitle="rain water mixing ratio [g/kg]",
            opacity=args.opacity2,
            xscalarbar=0.9,
            yscalarbar=0.1,
            logscale=args.logscale2,
            contour_value=args.contour2,
        )
        Hide(tempxmf, renderView1)
        hide_scalar_bar()
        time(renderView1, 0.15, 0.9)
        layout()
        current_camera_placement(args.camera)
        save_middle_frame(args.middleframename)
        SaveAnimation(args.animationname, renderView1, FrameRate = args.framerate)
    if args.real is True:
        display_settings(*args.boxcolor, *args.boxsize)
        background(args.photopath, photobackgroung=True)
        threshold1 = create_treshold(
            "Threshold1", tempxmf, "UnstructuredGridRepresentation", args.lowerthreshold
        )
        threshold2 = create_treshold(
            "Threshold2",
            tempxmf,
            "UnstructuredGridRepresentation",
            args.lowerthreshold,
            scalars=True,
        )
        calculator(
            registrationname="Calculator1",
            input=threshold1,
            function="rc*1000",
            resultarrayname="rc [g/kg]",
            colormap=args.colormap1,
            scalarbartitle="cloud water mixing ratio [g/kg]",
            opacity = args.opacity1,
            xscalarbar=0.9,
            yscalarbar=0.6,
            logscale=args.logscale1,
            contour_value=True,
            move_data=True,
        )
        calculator(
            registrationname="Calculator2",
            input=threshold2,
            function="rr*1000",
            resultarrayname="rr [g/kg]",
            colormap=args.colormap2,
            scalarbartitle="rain water mixing ratio [g/kg]",
            opacity=args.opacity2,
            xscalarbar=0.9,
            yscalarbar=0.1,
            logscale=args.logscale2,
            contour_value=True,
            move_data=True,
        )
        Hide(tempxmf, renderView1)
        hide_scalar_bar()
        if args.axisinrealistic is True:
            extract_bottom_xy()
        layout()
        current_camera_placement("Zoom")
        save_middle_frame(args.middleframename)
        SaveAnimation(args.animationname, renderView1, FrameRate=args.framerate)
