new Vue({
    el: '#app',
    vuetify: new Vuetify(),
    data: () => ({
        valid: true,
        name: '',
        nameRules: [
            v => !!v || 'Sequence is required',
            v => (v && v.length <= 20000) || 'Sequence must be less than 20,000 characters',
        ],
        email: '',
        emailRules: [
            v => !!v || 'E-mail is required',
            v => /.+@.+\..+/.test(v) || 'E-mail must be valid',
        ],
        select: null,
        items: [
            'Item 1',
            'Item 2',
            
            'Item 3',
            'Item 4',
        ],
        checkbox: false,
    }),
  
    methods: {
        validate () {
            if (this.$refs.form.validate()) {
                this.snackbar = true
            }
        },
        reset () {
            this.$refs.form.reset()
        },
        resetValidation () {
            this.$refs.form.resetValidation()
        },
    },
  })