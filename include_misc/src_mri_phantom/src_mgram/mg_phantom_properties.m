function phantom = mg_phantom_properties(phantom, varargin)

phantom.p1D.ind   = phantom.ind_map(phantom.mask2D);
phantom.p2D.ind   = phantom.ind_map;
phantom.p1D.pos   = phantom.pos_map(phantom.mask2D);
phantom.p2D.pos   = phantom.pos_map;
phantom.p1D.mask  = phantom.mask1D;
phantom.p2D.mask  = phantom.mask2D;
phantom.p1D.cmaps = phantom.cmaps(:,phantom.mask1D).';
phantom.p2D.cmaps = phantom.cmaps;

for j=1:numel(varargin)
    prop = varargin{j};
    if numel(prop)~=phantom.n_phant+1
        error(['wrong number in: ' inputname(j+1)]);
    end
    prop2D = zeros(phantom.Nxy, phantom.Nxy);
    for p = 1 : phantom.n_phant+1        
        prop2D(phantom.ind_map==p) = prop(p);
    end
    prop1D = prop2D(phantom.mask2D);
    eval(['phantom.p.' inputname(j+1) '=prop;']);
    eval(['phantom.p2D.' inputname(j+1) '=prop2D;']);
    eval(['phantom.p1D.' inputname(j+1) '=prop1D;']);
    clear prop1D prop2D;
end

phantom = rmfield(phantom, {'cmaps'; 'ind_map'; 'pos_map'; 'mask2D'; 'mask1D'});

end

