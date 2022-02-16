function print3( fig,fname )
print(fig,[fname '.eps'],'-depsc','-r300')
print(fig,[fname '.pdf'],'-dpdf','-r300')
print(fig,[fname '.png'],'-dpng','-r300')
end

