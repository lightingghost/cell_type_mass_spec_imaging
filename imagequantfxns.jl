
function detect_nuclei(img, σscales, minimal_blob_amp, gaussianblur) # σscales =  2.0.^[3]
  original_img = copy(img)
  img = imfilter(img,Kernel.gaussian(gaussianblur)) # gaussianblur = 3

  # σscales =  2.0.^[0.5,0,1]

  blobs = blob_LoG(img, σscales, true)

  amps =[]
  for blob in blobs
    amp = blob.amplitude
    push!(amps, amp)
  end
  println("Total detected blobs:")
  println(length(blobs))

  ## this is important for keeping down spurious blobs
  kept_blobs= []
  for blob in blobs
    if blob.amplitude > minimal_blob_amp 
      push!(kept_blobs, blob)
    end
  end

 kept_blobs,mindistance=remove_overlapping(kept_blobs)
 colorimg = RGB.(original_img)
  colorimg = copy(colorimg)


  circle_img = draw_circle_at_blob(colorimg, kept_blobs)
  return length(kept_blobs), kept_blobs, circle_img,original_img 
end

function get_blob_chars(imgfilename, blobs, images_dict, imgdict)
  img = (images_dict[imgfilename])
  new_img = copy(img)

  imgs = Array{Any}(undef, length(imgdict))
  for i in 1:length(imgdict)
    imgs[i] = new_img[i,:,:][1]
  end

  blob_chars = []
  for (i,blob) in enumerate(blobs)
    if i % 10 == 0
      println(i)
    end
      radius = (blob.σ)*(sqrt(2))
      center = blob.location
      temparray = []
      for ra in [0.5 1 1.2]
        sums = []
        graymask = make_mask(blob, imgs[1], ra)
        mask = Array{Float64}(graymask)
        sum_mask = sum(mask)
        for (k,v) in imgdict
          push!(sums, sum(mask .* imgs[imgdict[k]]))
       end
        push!(temparray, [sum_mask sums...])
      end
      temparray = hcat(blob, temparray...)
      push!(blob_chars, temparray)
    end
  return blob_chars
end

function getotsus(dir, fullpath, conf)
    fillsize = conf["fillsize"]
    gaussianblur = conf["gaussianblur"]
    minimal_blob_amp = conf["minimal_blob_amp"]
    σscales = conf["σscales"] 
    rev_channeldict = conf["rev_channeldict"]
    channeldict = conf["channeldict"]
	arrayout = []
    images_dict = Dict()
	for tif in readdir(fullpath) 
		if occursin(r".tif$", tif)
			images_dict[tif] = Array{Any}(undef, 4) 
			images_dict[tif][4] = "wt"
			imgfilename = joinpath(fullpath, tif)
			img = FileIO.load(imgfilename)
			for i in 1:3
				if i == 1
                    chanimg = map(x -> x.b, img)
				elseif i == 2
                    chanimg = map(x -> x.g, img)
				elseif i == 3
                    chanimg = map(x -> x.r, img)
				end
				chanimg = Gray.(chanimg)
				images_dict[tif][i] = parent(padarray(chanimg, Fill(0, (50,50), (50,50)))) # load channels, padding each with 50 pixels on a side
			end

			nuclei_count, nuclei, circle_img, dapi_only_img = detect_nuclei(images_dict[tif][1], σscales, minimal_blob_amp, gaussianblur)

			blob_chars = get_blob_chars(tif, nuclei, images_dict, rev_channeldict) 

			blob_chars = vcat(blob_chars...)
			imgfilecol = fill(imgfilename, size(blob_chars[:,2]))
			blob_chars = hcat(imgfilecol, blob_chars)

			# getting stats for cell types
			ki67core = blob_chars[:,6]
			gfapfull = blob_chars[:,9]
			thresh_gfap = otsu_threshold(Gray.(gfapfull))
			push!(arrayout, thresh_gfap)
		end
	end
	return arrayout
end

function remove_overlapping(kept_blobs)
    # look for overlap between blobs
    mindistance = []
    i = 1
    while i <= size(kept_blobs)[1]
        blob= kept_blobs[i]
        s = length(kept_blobs)
        distmatrix = Array{Float64,1}(undef, s)
       for (ii,blob2) in enumerate(kept_blobs)
           distmatrix[ii] = dist(blob.location,blob2.location)
       end
       push!(mindistance, minimum(distmatrix[map(x->x>0, distmatrix)]))
       distmatrix = map(x -> x < (sqrt(2)*blob.σ), distmatrix)
       i += 1
       if sum(distmatrix) > 1
           kept_blobs = deleteat!(kept_blobs,i)
           i = 1
       end
    end
    return kept_blobs, mindistance
end

function dist(x::CartesianIndex,y::CartesianIndex)
    return sqrt((x.I[1] - y.I[1])^2 + (x.I[2] - y.I[2])^2)
end

function draw_circle_at_blob(img_orig, blobs)
  img = copy(img_orig)
  ind = axes(img)
  for blob in blobs
    radius = (blob.σ)*(sqrt(2))
    center = blob.location
    for a in 0:(π*0.05):2π 
      for i in 0:0.01 
        try
          x = Int(round(center[1] + (radius+i)*cos(a)))
          y = Int(round(center[2] + (radius+i)*sin(a)))
          img[x,y] = RGB(1.0,0,0)
      catch
        end
        try
          x = Int(round(center[1] + (radius*1.2+i)*cos(a)))
          y = Int(round(center[2] + (radius*1.2+i)*sin(a)))
          img[x,y] = RGB(0,1.0,0)
      catch
        end
        try
          x = Int(round(center[1] + (radius*0.5+i)*cos(a)))
          y = Int(round(center[2] + (radius*0.5+i)*sin(a)))
          img[x,y] = RGB(0,0,1.0)
      catch
        end
      end
    end
    end
  return img
end

function make_mask(blob, img, ratio) 
    blankcanvas = zeros(size(img))
    blankcanvas = Gray.(blankcanvas)
    radius = (blob.σ)*(sqrt(2))
    radius = ratio*radius
    center = blob.location
    circ = ImageDraw.CirclePointRadius(center, radius)
    draw!(blankcanvas, circ)
    return blankcanvas
end

function process_tifs(dir, fullpath,fullpathout, conf)
    fillsize = conf["fillsize"]
    gaussianblur = conf["gaussianblur"]
    minimal_blob_amp = conf["minimal_blob_amp"]
    σscales = conf["σscales"] 
    rev_channeldict = conf["rev_channeldict"]
    channeldict = conf["channeldict"]
	arrayout = []
    images_dict = Dict()
	for tif in readdir(fullpath) 
		if occursin(r".tif$", tif)
			images_dict[tif] = Array{Any}(undef, 4) 
			images_dict[tif][4] = "wt"
			imgfilename = joinpath(fullpath, tif)
			img = FileIO.load(imgfilename)
			for i in 1:3
				if i == 1
                    chanimg = map(x -> x.b, img)
				elseif i == 2
                    chanimg = map(x -> x.g, img)
				elseif i == 3
                    chanimg = map(x -> x.r, img)
				end
				chanimg = Gray.(chanimg)
				images_dict[tif][i] = parent(padarray(chanimg, Fill(0, (50,50), (50,50)))) 
			end

			nuclei_count, nuclei, circle_img, dapi_only_img = detect_nuclei(images_dict[tif][1], σscales, minimal_blob_amp, gaussianblur)

			blob_chars = get_blob_chars(tif, nuclei, images_dict, rev_channeldict) 

			blob_chars = vcat(blob_chars...)
			imgfilecol = fill(imgfilename, size(blob_chars[:,2]))
			blob_chars = hcat(imgfilecol, blob_chars)

			# getting stats for cell types
			ki67core = blob_chars[:,6]
			gfapfull = blob_chars[:,9]
			thresh_gfap = otsu_threshold(Gray.(gfapfull))
			blob_chars = hcat(blob_chars, map(x-> x > thresh_gfap, gfapfull))

			header = ["file" "blobstats" "core_mask" "core_dapi" "core_gfap" "core_ki67" "full_mask" "full_dapi" "full_gfap" "full_ki67" "cyto_mask" "cyto_dapi" "cyto_gfap" "cyto_ki67"  "gfap_pos"]
			dataout = vcat(header, blob_chars)
            writedlm(joinpath(fullpathout, "$tif.csv"), dataout, ',')

            gfap_imgcircled, dapi_imgcircled = id_blobs(blob_chars, fullpathout, images_dict, tif)
            fifour = [gfap_imgcircled; dapi_imgcircled]
            FileIO.save(joinpath(fullpathout,"validation_gfapbottom_posneg_$tif.jpg"),fifour)
			#summary stats
            percent_gfap, percent_dapialone = breakdown_celltypes(blob_chars)
            gfap_count= Int(round(percent_gfap*nuclei_count))
            nongfap_count = nuclei_count - gfap_count
			push!(arrayout, [tif percent_gfap percent_dapialone nuclei_count gfap_count nongfap_count])
		end
	end
	return arrayout
end

function breakdown_celltypes(blob_chars)
	gfap = 0
	dapialone = 0
	total_nuclei = size(blob_chars)[1]
	for i in 1:size(blob_chars)[1] 
		gfappos = blob_chars[i,15]
		if gfappos
			gfap += 1
        else
			dapialone +=1
		end
	end
	percent_gfap_alone = gfap/total_nuclei
	percent_dapialone = dapialone/total_nuclei
	return percent_gfap_alone, percent_dapialone
end

function id_blobs(blob_chars, fullpathout, images_dict, tif)
    img_base = colorview(RGB, Float64.(images_dict[tif][3]), Float64.(images_dict[tif][2]), Float64.(images_dict[tif][1]))

    gfapposBlobs =[]
    dapialone = []
	for i in 1:size(blob_chars)[1] 
		blob = blob_chars[i,2]
		radius = (blob.σ)*(sqrt(2))
		center = blob.location
		gfappos = blob_chars[i,15]

        if (gfappos) 
            push!(gfapposBlobs, i)
        else
            push!(dapialone,i)
        end
    end
    gfap_imgcircled = copy(img_base)
    dapi_imgcircled = copy(img_base)
    if length(gfapposBlobs) >0
        for i in gfapposBlobs
            bloboi = blob_chars[i,2]
            gfap_imgcircled = drawcirclesimple(gfap_imgcircled, bloboi)
        end
    end
   if length(dapialone) >0
        for i in dapialone
            bloboi = blob_chars[i,2]
            dapi_imgcircled = drawcirclesimple(dapi_imgcircled,bloboi)
        end
    end
    return gfap_imgcircled, dapi_imgcircled 
end

function drawcirclesimple(img,blob)
		radius = (blob.σ)*(sqrt(2))
		center = blob.location
		for a in 0:(π*0.05):2π 
        for i in 0:0.01 
            try
                x = Int(round(center[1] + (radius+i)*cos(a)))
                y = Int(round(center[2] + (radius+i)*sin(a)))
                img[x,y] = RGB(1.0,1.0,1.0)
            catch
            end
        end
    end
    return img
end
