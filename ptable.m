function varargout = ptable(varargin)
%PTABLE(Z,A) Determines mass of given atom using semi-empirical formula
%   Z = #protons, A = #nucleons
%
%   0 args:                   plots mass, binding energy v. atomic number
%   1 args, Z:                returns atomic mass & struct of atomic information
%   2 args, Z & A:            returns the mass of isotope & whether it is stable
%   2 args, Z & field name:   returns atomic mass & value of field name

ptabledata = readtable('periodictabledata1.csv');

mp = 1.6726219*10^-27; %proton mass, kg
mn = 1.674929*10^-27; %neutron mass, kg

switch nargin
    %% When no arguments: 1) SEmass, atomic weight v. atomic number 2) binding E per A v. atomic number
    case 0 
        %Create arrays of Z, A, atomic weights, SE masses, and binding energies
        zvalues = ptabledata.AtomicNumber;
        aweights = ptabledata.AtomicWeight.*(931493614.838475/(10^6));  %convert amu to MeV
        avalues = round(ptabledata.AtomicWeight);
        [semvalues, evalues] = semass(avalues, zvalues);

        %Determine binding energy per nucleon
        eper = evalues./avalues;

        %Plot SEM, atomic weight v. Z
        clf
        figure(1), hold on;
            plot(zvalues, semvalues, 'b--')
            plot(zvalues, aweights, 'k--')
            legend('Semi-empirical mass (MeV)', 'Atomic weight (MeV)')
            xlabel('Atomic Number')
            ylabel('Semi-empirical mass (MeV)') %something still wrong with this plot
            title('Semi-empirical Mass v. Atomic Number')
            savefig('sem.fig')

        %Plot binding energy per nucleon v. Z
        figure(2)
            plot(zvalues, eper, 'b-')
            xlabel('Atomic Number')
            ylabel('Binding energy per nucleon (MeV)')
            title('Binding Energy Per Nucleon v. Atomic Number')
            savefig('eperA.fig')

    %% when 1 num arg (z), return: 1) SE mass of element, 2) structure w/ all fields for element
    case 1
        Z = varargin{1};

        aweight = ptabledata.AtomicWeight(Z);
        A = round(aweight);
        profile.labels = ptabledata.Properties.VariableNames;
        profile.data = ptabledata(Z,:);

        %calculate semi-empirical mass
        [se_mass, ~] = semass(A,Z);
        
        varargout{1} = se_mass;
        varargout{2} = profile;

    
        
    %% 2 cases for when 2 input values:
    case 2
       % when 2 num args (z and # nucleons): return 1)mass of isotope, 2)stability
        if isnumeric(varargin{1}) && isnumeric(varargin{2})
            Z = varargin{1};
            A = varargin{2};
            
            [~, eb] = semass(A,Z);

            %Calculate mass of isotope
            isomass = (mp*Z + mn*(A-Z))*6.022*10^26; %converted from kg to amu

            %Determine stability; is this true?
            if Z>82 || eb<0 %1~=(A-Z)/Z
                stable = 'False';
            else
                stable = 'True';
            end
            
            varargout{1} = isomass;
            varargout{2} = stable;
            
       %when 1 num arg (z) and 1 char arg (field name in csv file): return 1)SEmass, 2)value of given field
        elseif isnumeric(varargin{1}) && ischar(varargin{2})
            Z = varargin{1};
            field = varargin{2};

            %grab value contained in given column
            fieldvalue = ptabledata(Z,field);
            aweight = ptabledata.AtomicWeight(Z);
            A = round(aweight);

            %calculate SE mass of element
            [se_mass, ~] = semass(A,Z);
            
            varargout{1} = se_mass;
            varargout{2} = fieldvalue;
        end   
end
end
