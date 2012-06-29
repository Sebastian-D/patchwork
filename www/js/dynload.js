$(document).ready(function(){  

    //References  
    var sections = $("#menu li");  
    var loading = $("#loading");  
    var content = $("#content");  

    //hide loading bar  
function hideLoading(){  
    loading.fadeTo(1000, 0);  
    content.slideDown();  
}; 

//Manage click events  
sections.click(function(){  
    //show the loading bar  
    showLoading();  
    //load selected section  
    switch(this.id){  
        case "install":  
            content.slideUp();  
            content.load("test.php #section_installation", hideLoading);  
            break;  
        case "execution":  
            content.slideUp();  
            content.load("test.php #section_execution", hideLoading);  
            break;  
        case "results":  
            content.slideUp();  
            content.load("test.php #section_results", hideLoading);  
            break;  
        default:  
            //hide loading bar if there is no selected section  
            hideLoading();  
            break;  
    }  
});  

    });  