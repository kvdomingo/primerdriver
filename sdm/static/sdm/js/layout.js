document.addEventListener('DOMContentLoaded', () => {
    const cl = new cloudinary.Cloudinary({
        cloud_name: 'kdphotography-assets',
        secure: true,
    });
    cl.responsive();

    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    });
});
