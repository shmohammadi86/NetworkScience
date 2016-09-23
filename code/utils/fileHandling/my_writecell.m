function [ res ] = my_writecell( path, input_cell, dlm)
    function writeElement(element)
        if(isempty(element))
             fprintf(fd, '-');
        else            
            if(ischar(element))
                fprintf(fd, '%s', element);
            elseif(isnumeric(element))
                if(floor(element)==element)
                    fprintf(fd, '%d', element);                    
                else
                    temp = shift*element;
                    if(floor(temp) == temp)
                        fprintf(fd, '%f', element);
                    else
                        fprintf(fd, '%e', element);
                    end
                end
            else
                error('my_writecell:: Unknown datatype for element (%d, %d) in the cell', i, j);
            end                       
        end
    end

    if(nargin < 2)
        error('my_writecell:: Invalid number of input argument (%d)', nargin);
    elseif(nargin == 2)
        dlm = '\t';
    end
    
    decimal_points = 5;
    shift = 10^decimal_points;
    
    if(isempty(input_cell))
        res = -1;
        return;
    end
    
    fd = fopen(path, 'w');
    if(fd == -1)
        error('my_writecell:: Cannot open file %s', path);
    end
    
    m = size(input_cell, 1);
    n = size(input_cell, 2);
    for i = 1:m
        writeElement(input_cell{i, 1});
        for j = 2:n
            fprintf(fd, dlm);
            writeElement(input_cell{i, j});           
        end
        fprintf(fd, '\n');
    end
    
    fclose(fd);
    res = 0;    
end