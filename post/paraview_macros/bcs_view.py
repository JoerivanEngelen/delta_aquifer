#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'NetCDF Reader'\bcs.nc
path=r"c:\Users\engelen\test_imodpython\synth_delta_test\SD_i239\input\data\bcs.nc"
bcsnc = NetCDFReader(FileName=[path])
bcsnc.Dimensions = '(z, y, x)'

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1009, 857]

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
#transform2.Transform.Translate = [0.0, 150000.0, 0.0]
transform2.Transform.Translate = [100.0, -1000.0, 1.0]
transform2.Transform.Scale = [1.0, 1.0, 100.0]

# show data in view
transform2Display = Show(transform2, renderView1)
# trace defaults for the display properties.
transform2Display.Representation = 'Outline'
transform2Display.ColorArrayName = [None, '']
transform2Display.OSPRayScaleArray = 'river_stage'
transform2Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform2Display.SelectOrientationVectors = 'None'
transform2Display.ScaleFactor = 19900.0
transform2Display.SelectScaleArray = 'None'
transform2Display.GlyphType = 'Arrow'
transform2Display.GlyphTableIndexArray = 'None'
transform2Display.DataAxesGrid = 'GridAxesRepresentation'
transform2Display.PolarAxes = 'PolarAxesRepresentation'
transform2Display.ScalarOpacityUnitDistance = 1845.2831649390853

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 100.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(programmableFilter2, renderView1)

# update the view to ensure updated data information
renderView1.Update()

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

# set scalar coloring
ColorBy(threshold2Display, ('CELLS', 'river_stage'))

# rescale color and/or opacity maps used to include current data range
threshold2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(threshold2, renderView1)

ColorBy(threshold3Display, ('CELLS', 'river'))

riverLUT = GetColorTransferFunction('river')
riverLUT.RGBPoints = [0.0, 0.9, 0.9, 0.9, 2.0, 0.9, 0.9, 0.9]
riverLUT.ScalarRangeInitialized = 1.0
riverLUT.EnableOpacityMapping = 1
riverLUT.ScalarOpacityFunction = CreatePiecewiseFunction(Points=[0.0, 0.5, 0.5, 0.0, 100.0, 0.5, 0.5, 0.0])

threshold3 = Threshold(Input=transform2)
threshold3.Scalars = ['CELLS', 'sea']
threshold3.ThresholdRange = [0.17, 1.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1009, 857]

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
lithLUT.ScalarOpacityFunction = CreatePiecewiseFunction(Points=[0.0, 0.5, 0.5, 0.0, 100.0, 0.5, 0.5, 0.0])


#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [442684.2037462758, -406480.74239496514, 224788.59000526415]
renderView1.CameraFocalPoint = [99999.99999999993, -3.7612197757863344e-11, -440.335021972661]
renderView1.CameraViewUp = [-0.3056007632794434, 0.2512504923945124, 0.9184124147432551]
renderView1.CameraParallelScale = 149441.3021946477



#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).