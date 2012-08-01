
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="http://r-forge.r-project.org/themes/rforge//images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p> This project has been launched to fit mixture of expert models and simultaneously select important variables through regularization. The related statistical articles are</p>
<p>Khalili, A. and Chen, J. (2007). Variable Selection in Finite Mixture of Regression Models.
        Journal of the American Statistical Association, 102, 1025-1038. </p>
<p>Khalili, A. (2010). New Estimation and Feature Selection Methods in Mixture-of-Experts Models.
        The Canadian Journal of Statistics, 38, 519-539. </p> 
<p>The current developers are <a href="http://vahid.probstat.ca/"><strong>Vahid Partovi Nia</strong></a>, 
<a href="http://www.math.mcgill.ca/khalili/"><strong>Abbas Khalili</strong></a> and <a href="delphine.html"><strong>Delphine Savignard</strong></a>. </p>
<div style="text-align: center"> 
<img src="expert1.png" width="589" height="558" alt="" />
</div>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
