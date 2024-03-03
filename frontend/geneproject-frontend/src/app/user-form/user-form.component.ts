import { Component } from '@angular/core';
import {FormBuilder, Validators} from '@angular/forms';
import {DropdownModule} from 'primeng/dropdown';
import { MultiSelectModule } from 'primeng/multiselect';
@Component({
  selector: 'app-user-form',
  templateUrl: './user-form.component.html',
  styleUrl: './user-form.component.css'
})

export class UserFormComponent {
  ethnicityOptions = [
    {name:'Asian'}, 
    {name:'African'},
    {name:'Caucasian'},
    {name: 'Hispanic/Latino'}
  ]
  
  userForm = this.fb.group({ 
    firstName: ['', Validators.required],
    lastName: ['', Validators.required],
    ethnicity: ['', Validators.required],
    height : ['', Validators.required],
    weight : ['', Validators.required],
    budget : [''],
  });
  constructor(public fb: FormBuilder){
  
  }
  submit(){

  }
}
