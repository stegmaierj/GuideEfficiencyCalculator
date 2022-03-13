function [] = export_png_without_border(filename, figureHandle, width, height)
    set(figureHandle, 'PaperPositionMode', 'Manual')
    set(figureHandle, 'PaperUnits', 'Points')
    set(figureHandle, 'PaperPosition', [0,0,width, height])
    print(figureHandle, filename, '-dpng')
end