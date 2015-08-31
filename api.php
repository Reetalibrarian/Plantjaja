<?php
session_start();

ini_set("upload_max_filesize", "20M");

$phyllotaxy = empty($_REQUEST['phyllotaxy'])?1:$_REQUEST['phyllotaxy'];

$session_id = md5(microtime(true));

$target_dir = "/var/www/html/plantjaja/uploads/$session_id/";
if (!is_dir($target_dir)) {
	mkdir($target_dir, 0777, true);
}

$target_file = $target_dir . basename($_FILES["fileToUpload"]["name"]);
$type = $_FILES["fileToUpload"]["type"];
$uploadOk = 1;
$imageFileType = pathinfo($target_file,PATHINFO_EXTENSION);
// Check if image file is a actual image or fake image
if(isset($_POST["submit"])) {
	$check = getimagesize($_FILES["fileToUpload"]["tmp_name"]);
	if($check !== false) {
		$uploadOk = 1;
	}
	else {
		$uploadOk = 0;
	}
}

if ($imageFileType != "jpg" && $imageFileType != "jpeg" && $imageFileType != "JPG" && $imageFileType != "JPEG" && $type != 'image/jpeg' && $_FILES["fileToUpload"]["name"] != 'blob') {
	$uploadOk = 0;
}
else if ($type == 'image/jpeg' && $_FILES["fileToUpload"]["name"] == 'blob') {
	$uploadOk = 1;
}

if ($uploadOk) {
	if (move_uploaded_file($_FILES["fileToUpload"]["tmp_name"], $target_file)) {

		/*
		list($width, $height) = getimagesize($target_file);
		$thumb = imagecreatetruecolor(720, 960);
		$source = imagecreatefromjpeg($target_file);
		imagecopyresized($thumb, $source, 0, 0, 0, 0, 720, 960, $width, $height);
		imagejpeg($thumb, $target_file);
		//*/

		// run matlab program
		// echo "majaja";
		$p_file_tpl = "Plant_classifier_ouline_folder_phyllotaxy_2.m.tpl";
		$l_file = "Leaf_Feature_extract.m";
		$conf_tpl = "conf.m.tpl";
		#$workspace_base = "/usr/local/MATLAB/R2012a/mmder/plantjaja/workspace/";
		$workspace_base = "/var/www/html/plantjaja/workspace/";
		$exec_dir = $workspace_base . $session_id . "/";
		if (!is_dir($exec_dir)) {
			mkdir($exec_dir, 0777, true);
		}

		$conf_content = str_replace("{exec_dir}", $exec_dir, file_get_contents($conf_tpl));
		$conf_content = str_replace("{upload_dir}", $target_dir, $conf_content);
		$conf_content = str_replace("{filename}", basename($target_file), $conf_content);
		$conf_content = str_replace("{phyllotaxy}", $phyllotaxy, $conf_content);

		$p_content = file_get_contents($p_file_tpl);

		file_put_contents($exec_dir . basename($p_file_tpl, '.tpl'), $conf_content . "\n" . $p_content);
		copy($l_file, $exec_dir . basename($l_file));
		exec("cp /var/www/html/plantjaja/dep/* " . $exec_dir);

		exec('export MATLAB_PREFDIR=/var/www/html/plantjaja/.matlab && /usr/local/MATLAB/R2012a/bin/matlab -nodesktop -nosplash -nodisplay -r "run '. $exec_dir . basename($p_file_tpl, '.tpl') . '; quit;"');

		$o = explode("  ", trim(file_get_contents($target_dir . "Sample_0001.txt"), " \r\n"));
		/*
		echo "<xmp>";
		var_dump($o);
		echo "</xmp>";
		//*/
		$o_translated = array_map("idToName", $o);
		echo json_encode($o_translated);

	}
	else {
		echo "move file failed\n";
	}
}
else {
	echo "<xmp>";
	var_dump($_FILES);
	echo "</xmp>";
}

session_destroy();

function idToName ($id) {
	$id = trim($id);
	$f = file("SpeciesNum.txt");
	$map = array();
	foreach ($f as $line) {
		$line = trim($line, " \r\n");
		$parts = explode("\t", $line);
		$map[$parts[0]] = $parts[1];
	}
	return $map[$id];
}



?>
