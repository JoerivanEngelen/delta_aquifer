#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import glob
import os, sys
import math

# create a new 'NetCDF Reader'
fol=sys.argv[1]
nc_paths = glob.glob(os.path.join(fol, "*[0-9][0-9][0-9].nc"))
nc_paths.sort()

netCDFReader1 = NetCDFReader(FileName=nc_paths)
netCDFReader1.Dimensions = '(z, y, x)'

bc_path=sys.argv[2]
bcsnc = NetCDFReader(FileName=[bc_path])
bcsnc.Dimensions = '(z, y, x)'

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
calculator1.ResultArrayName = 'v'
calculator1.Function = '(-1*vx)*iHat+(1*vy)*jHat+(1*vz)*kHat'

# hide data in view
Hide(netCDFReader1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=calculator1)

# Properties modified on programmableFilter1
programmableFilter1.Script='inp = self.GetInput()\ndims = inp.GetDimensions()\next = inp.GetExtent()\n\noup = self.GetOutput()\noup.SetDimensions(dims[0]+1, dims[1]+1, dims[2]+1)\noup.SetExtent(ext[0], ext[1]+1, ext[2], ext[3]+1, ext[4], ext[5]+1)\n\nN=inp.GetPointData().GetNumberOfArrays()\n\nfor i in range(N):\n    data = inp.GetPointData().GetAbstractArray(i)\n    oup.GetCellData().AddArray(data)\n\n'
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# hide data in view
Hide(calculator1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

#Get z
z_min, z_max = calculator1.GetDataInformation().GetBounds()[4:]
dz = z_max - z_min

zscale = 4000. / math.sqrt(dz)

# create a new 'Transform'
transform1 = Transform(Input=programmableFilter1)
transform1.Transform = 'Transform'

# Properties modified on transform1.Transform
transform1.Transform.Scale = [1.0, 1.0, zscale]

# show data in view
transform1Display = Show(transform1, renderView1)
# trace defaults for the display properties.
transform1Display.Representation = 'Outline'
transform1Display.ColorArrayName = [None, '']
transform1Display.OSPRayScaleArray = 'conc1'
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
threshold1.Scalars = ['CELLS', 'conc1']
threshold1.ThresholdRange = [-1.0, 37]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = [None, '']
threshold1Display.OSPRayScaleArray = 'conc1'
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
ColorBy(threshold1Display, ('CELLS', 'conc1'))

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=threshold1)

last_t = animationScene1.TimeKeeper.TimestepValues[-1]
scale = 1./365.

annotateTimeFilter1.Format = 'Time: %.0f BP'
annotateTimeFilter1.Shift = last_t*scale
annotateTimeFilter1.Scale = -scale

annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

# get color transfer function/color map for 'conc1'
conc1LUT = GetColorTransferFunction('conc1')
conc1LUT.RGBPoints = [-0.018854578956961632, 0.171875, 0.48046875, 0.7109375, 0.9965649232740552, 0.171875, 0.48046875, 0.7109375, 0.9965649232740552, 0.5703125, 0.7734375, 0.87109375, 5.058242932198123, 0.5703125, 0.7734375, 0.87109375, 5.058242932198123, 0.8671875, 0.9375, 0.8125, 10.135340443353208, 0.8671875, 0.9375, 0.8125, 10.135340443353208, 0.99609375, 0.87109375, 0.6015625, 15.212437954508292, 0.99609375, 0.87109375, 0.6015625, 15.212437954508292, 0.9609375, 0.5625, 0.32421875, 20.289535465663377, 0.9609375, 0.5625, 0.32421875, 20.289535465663377, 0.83984375, 0.09765625, 0.109375, 25.36663297681846, 0.83984375, 0.09765625, 0.109375, 25.36663297681846, 0.62890625, 0.09765625, 0.109375, 30.443730487973546, 0.62890625, 0.09765625, 0.109375, 30.443730487973546, 0.421875, 0.09765625, 0.109375, 36.02853775024414, 0.421875, 0.09765625, 0.109375]
conc1LUT.ScalarRangeInitialized = 1.0

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Programmable Filter'
programmableFilter2 = ProgrammableFilter(Input=bcsnc)
programmableFilter2.Script = ''
programmableFilter2.RequestInformationScript = ''
programmableFilter2.RequestUpdateExtentScript = ''
programmableFilter2.PythonPath = ''

# set active source
SetActiveSource(programmableFilter2)

# Properties modified on programmableFilter2
programmableFilter2.OutputDataSetType = 'vtkImageData'
programmableFilter2.Script = 'dims = inputs[0].GetDimensions()\next = inputs[0].GetExtent()\noutput.SetDimensions(dims[0]+1, dims[1]+1, dims[2]+1)\noutput.SetExtent(ext[0], ext[1]+1, ext[2], ext[3]+1, ext[4], ext[5]+1)\ninputPd = inputs[0].PointData\noutputCd = output.CellData\nfor array in inputPd:\n   print(type(array))\n   outputCd.append(array, array.GetName())'
programmableFilter2.RequestInformationScript = ''
programmableFilter2.RequestUpdateExtentScript = ''
programmableFilter2.PythonPath = ''

# hide data in view
Hide(bcsnc, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Transform'
transform2 = Transform(Input=programmableFilter2)
transform2.Transform = 'Transform'

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=transform2.Transform)

# Properties modified on transform2.Transform
transform2.Transform.Translate = [100.0, -1000.0, 1.0]
transform2.Transform.Scale = [1.0, 1.0, zscale]

# hide data in view
Hide(programmableFilter2, renderView1)

# create a new 'Threshold'
threshold2 = Threshold(Input=transform2)
threshold2.Scalars = ['CELLS', 'river']
threshold2.ThresholdRange = [0.1, 1]

# show data in view
threshold2Display = Show(threshold2, renderView1)
# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'
threshold2Display.ColorArrayName = [None, '']
threshold2Display.OSPRayScaleArray = 'river_stage'
threshold2Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold2Display.SelectOrientationVectors = 'None'
threshold2Display.ScaleFactor = 17200.0
threshold2Display.SelectScaleArray = 'None'
threshold2Display.GlyphType = 'Arrow'
threshold2Display.GlyphTableIndexArray = 'None'
threshold2Display.DataAxesGrid = 'GridAxesRepresentation'
threshold2Display.PolarAxes = 'PolarAxesRepresentation'
threshold2Display.ScalarOpacityUnitDistance = 2957.441091203209
threshold2Display.GaussianRadius = 8600.0
threshold2Display.SetScaleArray = [None, '']
threshold2Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold2Display.OpacityArray = [None, '']
threshold2Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

ColorBy(threshold2Display, ('CELLS', 'river'))

riverLUT = GetColorTransferFunction('river')
riverLUT.RGBPoints = [0.0, 0.9, 0.9, 0.9, 2.0, 0.9, 0.9, 0.9]
riverLUT.ScalarRangeInitialized = 1.0
riverLUT.EnableOpacityMapping = 1
riverLUT.ScalarOpacityFunction = CreatePiecewiseFunction(Points=[0.0, 0.5, 0.5, 0.0, 100.0, 0.5, 0.5, 0.0])

threshold3 = Threshold(Input=transform2)
threshold3.Scalars = ['CELLS', 'sea']
threshold3.ThresholdRange = [0.17, 1.0]

# show data in view
threshold3Display = Show(threshold3, renderView1)
# trace defaults for the display properties.
threshold3Display.Representation = 'Surface'
threshold3Display.ColorArrayName = [None, '']
threshold3Display.OSPRayScaleArray = 'river_stage'
threshold3Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold3Display.SelectOrientationVectors = 'None'
threshold3Display.ScaleFactor = 17700.0
threshold3Display.SelectScaleArray = 'None'
threshold3Display.GlyphType = 'Arrow'
threshold3Display.GlyphTableIndexArray = 'None'
threshold3Display.DataAxesGrid = 'GridAxesRepresentation'
threshold3Display.PolarAxes = 'PolarAxesRepresentation'
threshold3Display.ScalarOpacityUnitDistance = 8889.66660598137
threshold3Display.GaussianRadius = 8850.0
threshold3Display.SetScaleArray = [None, '']
threshold3Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold3Display.OpacityArray = [None, '']
threshold3Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold3Display.OpacityTransferFunction.Points = [0.0, 0.5, 0.5, 0.0, 100.0, 0.5, 0.5, 0.0]

ColorBy(threshold3Display, ('CELLS', 'sea'))

seaLUT = GetColorTransferFunction('sea')
seaLUT.RGBPoints = [0.0, 0.06, 0.06, 0.06, 2.0, 0.06, 0.06, 0.06]
seaLUT.ScalarRangeInitialized = 1.0
seaLUT.EnableOpacityMapping = 1
seaLUT.ScalarOpacityFunction = CreatePiecewiseFunction(Points=[0.0, 0.5, 0.5, 0.0, 100.0, 0.5, 0.5, 0.0])

threshold4 = Threshold(Input=transform2)
threshold4.Scalars = ['CELLS', 'lith']
threshold4.ThresholdRange = [1.5, 100.0]

# show data in view
threshold4Display = Show(threshold4, renderView1)
# trace defaults for the display properties.
threshold4Display.Representation = 'Surface'
threshold4Display.ColorArrayName = [None, '']
threshold4Display.OSPRayScaleArray = 'river_stage'
threshold4Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold4Display.SelectOrientationVectors = 'None'
threshold4Display.ScaleFactor = 17700.0
threshold4Display.SelectScaleArray = 'None'
threshold4Display.GlyphType = 'Arrow'
threshold4Display.GlyphTableIndexArray = 'None'
threshold4Display.DataAxesGrid = 'GridAxesRepresentation'
threshold4Display.PolarAxes = 'PolarAxesRepresentation'
threshold4Display.ScalarOpacityUnitDistance = 8889.66660598137
threshold4Display.GaussianRadius = 8850.0
threshold4Display.SetScaleArray = [None, '']
threshold4Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold4Display.OpacityArray = [None, '']
threshold4Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold4Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold4Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold4Display.OpacityTransferFunction.Points = [0.0, 0.5, 0.5, 0.0, 100.0, 0.5, 0.5, 0.0]

ColorBy(threshold4Display, ('CELLS', 'lith'))

lithLUT = GetColorTransferFunction('lith')
lithLUT.RGBPoints = [0.0, 0.6, 0.6, 0.6, 6.0, 0.6, 0.6, 0.6]
lithLUT.ScalarRangeInitialized = 1.0
lithLUT.EnableOpacityMapping = 1
lithLUT.ScalarOpacityFunction = CreatePiecewiseFunction(Points=[0.0, 0.2, 0.2, 0.0, 100.0, 0.2, 0.2, 0.0])

# update the view to ensure updated data information
renderView1.Update()

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
