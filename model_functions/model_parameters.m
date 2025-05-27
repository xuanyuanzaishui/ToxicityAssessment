function [p,c] = model_parameters(celltype)
p.celltype = celltype;

AF=0;
% channel conductances
c.G.GNa=23*(1-0.1*AF);   %(130)
c.G.GNaL=0.0025*AF; %(244)
c.G.Gto=0.165*(1.0-0.1*AF); %(341)
c.G.GKr=0.035; %(273)
c.G.GKs=0.0035*(1+0.5*AF); %(291)
c.G.GKur=0.045*(1.0-0.22*AF); %(363)
c.G.GK1=0.0525*(1+0.26*AF); %(384)
c.G.GCaL=(1-0.33*AF);%(147)
c.G.Pnak=1.26; %(132)
c.G.Gncx=3.15*(1+0.09*AF); %(152)
c.G.GpCa=0.0471;%(161)
c.G.SERCA = 1 ;%(456)
c.G.Leak = 1*(1.0+0.25*AF);%(457)
c.G.Rel = 1;
end
