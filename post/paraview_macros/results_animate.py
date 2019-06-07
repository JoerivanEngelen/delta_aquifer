#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import glob
import os, sys

# create a new 'NetCDF Reader'
fol=sys.argv[1]
nc_paths = glob.glob(os.path.join(fol, "*[0-9][0-9][0-9].nc"))
print(nc_paths)

netCDFReader1 = NetCDFReader(FileName=nc_paths)
netCDFReader1.Dimensions = '(z, y, x)'

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1003, 853]
# Properties modified on renderView1
renderView1.Background = [0.43137254901960786, 0.43137254901960786, 0.43137254901960786]


# create a new 'Calculator'
calculator1 = Calculator(Input=netCDFReader1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'v'
calculator1.Function = '(-1*vx)*iHat+(1*vy)*jHat+(1*vz)*kHat'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1418, 857]

# show data in view
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.Representation = 'Outline'
calculator1Display.ColorArrayName = [None, '']
calculator1Display.OSPRayScaleArray = 'conc'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'v'
calculator1Display.ScaleFactor = 19800.0
calculator1Display.SelectScaleArray = 'None'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'None'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.GaussianRadius = 9900.0
calculator1Display.SetScaleArray = ['POINTS', 'conc']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'conc']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(netCDFReader1, renderView1)

# find source
programmableFilter2 = FindSource('ProgrammableFilter2')

# find source
calculator2 = FindSource('Calculator2')

# find source
transform2 = FindSource('Transform2')

# find source
threshold2 = FindSource('Threshold2')

# find source
glyph2 = FindSource('Glyph2')

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=calculator1)
programmableFilter1.Script = ''
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# Properties modified on programmableFilter1
programmableFilter1.Script = 'dims = inputs[0].GetDimensions()\next = inputs[0].GetExtent()\noutput.SetDimensions(dims[0]+1, dims[1]+1, dims[2]+1)\noutput.SetExtent(ext[0], ext[1]+1, ext[2], ext[3]+1, ext[4], ext[5]+1)\ninputPd = inputs[0].PointData\noutputCd = output.CellData\nfor array in inputPd:\n   outputCd.append(array, array.GetName())\n'
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# show data in view
programmableFilter1Display = Show(programmableFilter1, renderView1)
# trace defaults for the display properties.
programmableFilter1Display.Representation = 'Outline'
programmableFilter1Display.ColorArrayName = [None, '']
programmableFilter1Display.OSPRayScaleArray = 'conc'
programmableFilter1Display.OSPRayScaleFunction = 'PiecewiseFunction'
programmableFilter1Display.SelectOrientationVectors = 'None'
programmableFilter1Display.ScaleFactor = 19800.0
programmableFilter1Display.SelectScaleArray = 'None'
programmableFilter1Display.GlyphType = 'Arrow'
programmableFilter1Display.GlyphTableIndexArray = 'None'
programmableFilter1Display.DataAxesGrid = 'GridAxesRepresentation'
programmableFilter1Display.PolarAxes = 'PolarAxesRepresentation'
programmableFilter1Display.GaussianRadius = 9900.0
programmableFilter1Display.SetScaleArray = [None, '']
programmableFilter1Display.ScaleTransferFunction = 'PiecewiseFunction'
programmableFilter1Display.OpacityArray = [None, '']
programmableFilter1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
programmableFilter1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
programmableFilter1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
programmableFilter1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(calculator1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Transform'
transform1 = Transform(Input=programmableFilter1)
transform1.Transform = 'Transform'

# Properties modified on transform1.Transform
transform1.Transform.Scale = [1.0, 1.0, 100.0]

# Properties modified on transform1.Transform
transform1.Transform.Scale = [1.0, 1.0, 100.0]

# show data in view
transform1Display = Show(transform1, renderView1)
# trace defaults for the display properties.
transform1Display.Representation = 'Outline'
transform1Display.ColorArrayName = [None, '']
transform1Display.OSPRayScaleArray = 'conc'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 19800.0
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityUnitDistance = 2225.3234597470287

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(programmableFilter1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Threshold'
threshold1 = Threshold(Input=transform1)
threshold1.Scalars = ['CELLS', 'conc']
threshold1.ThresholdRange = [-9999.0, 36.02853775024414]

# Properties modified on threshold1
threshold1.ThresholdRange = [-1.0, 36.02853775024414]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = [None, '']
threshold1Display.OSPRayScaleArray = 'conc'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 19500.0
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityUnitDistance = 3279.5243643393874
threshold1Display.GaussianRadius = 9750.0
threshold1Display.SetScaleArray = [None, '']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = [None, '']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(threshold1Display, ('CELLS', 'conc'))

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'conc'
concLUT = GetColorTransferFunction('conc')
concLUT.RGBPoints = [-0.018854578956961632, 0.171875, 0.48046875, 0.7109375, 0.9965649232740552, 0.171875, 0.48046875, 0.7109375, 0.9965649232740552, 0.5703125, 0.7734375, 0.87109375, 5.058242932198123, 0.5703125, 0.7734375, 0.87109375, 5.058242932198123, 0.8671875, 0.9375, 0.8125, 10.135340443353208, 0.8671875, 0.9375, 0.8125, 10.135340443353208, 0.99609375, 0.87109375, 0.6015625, 15.212437954508292, 0.99609375, 0.87109375, 0.6015625, 15.212437954508292, 0.9609375, 0.5625, 0.32421875, 20.289535465663377, 0.9609375, 0.5625, 0.32421875, 20.289535465663377, 0.83984375, 0.09765625, 0.109375, 25.36663297681846, 0.83984375, 0.09765625, 0.109375, 25.36663297681846, 0.62890625, 0.09765625, 0.109375, 30.443730487973546, 0.62890625, 0.09765625, 0.109375, 30.443730487973546, 0.421875, 0.09765625, 0.109375, 36.02853775024414, 0.421875, 0.09765625, 0.109375]
concLUT.ScalarRangeInitialized = 1.0

# create a new 'Glyph'
glyph1 = Glyph(Input=threshold1,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'None']
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleFactor = 19500.0
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.Vectors = ['CELLS', 'v']
glyph1.ScaleFactor = 19500.0
glyph1.GlyphMode = 'Every Nth Point'
glyph1.Stride = 500

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleArray = 'GlyphVector'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 20444.671289062502
glyph1Display.SelectScaleArray = 'GlyphVector'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'GlyphVector'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 10222.335644531251
glyph1Display.SetScaleArray = ['POINTS', 'conc']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'conc']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
#Hide(glyph2, renderView1)

# hide data in view
#Hide(threshold1, renderView1)

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'conc'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [438436.2873319544, -251656.11272326164, 134633.79718652414]
renderView1.CameraFocalPoint = [100500.00000000006, 38000.00000000013, -67477.78125000009]
renderView1.CameraViewUp = [-0.3474287370450981, 0.22912595429529117, 0.9092824477265138]
renderView1.CameraParallelScale = 126517.4334821998

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
SaveAnimation(os.path.join(fol, "screenshots", "results.png"), renderView1)